"""
create_users.py — Provision SOLANGE platform login accounts.

Creates a Supabase Auth user (email + password, email auto-confirmed) and a
matching public.users_profile row, so the person can sign in at
https://solange-platform.bio with their email and password.

Requires the Supabase SERVICE_ROLE key (admin privileges). This key is NEVER
committed — it is read from the environment at runtime. Get it from:
    Supabase dashboard → Project Settings → API → service_role (secret)

Usage
-----
    export SUPABASE_URL="https://lzzuxtnubznrkxwxjaab.supabase.co"
    export SUPABASE_SERVICE_KEY="<service_role secret key>"

    # Password auto-generated and printed once (recommended):
    python backend/create_users.py --email molly@cgu.edu --name "Molly Cohen"

    # Or set a specific password:
    python backend/create_users.py --email itamar.shabtai@cgu.edu \
        --name "Itamar Shabtai" --password "SomeStrongPass123!"

    # Make someone an admin (sees all runs + audit log):
    python backend/create_users.py --email me@cgu.edu --name "Doron Cohen" --role admin

Each run provisions one user. Re-running for an existing email is safe: the
script detects the existing account and only ensures the profile row exists.
"""

import argparse
import os
import secrets
import string
import sys

try:
    from supabase import create_client
except ImportError:
    sys.exit("supabase package not installed. Run: pip install supabase==2.4.2")


def generate_password(length: int = 16) -> str:
    """Strong random password: letters, digits, and a few safe symbols."""
    alphabet = string.ascii_letters + string.digits + "!@#$%*-_"
    # Guarantee at least one of each class so it passes any policy.
    while True:
        pw = "".join(secrets.choice(alphabet) for _ in range(length))
        if (any(c.islower() for c in pw) and any(c.isupper() for c in pw)
                and any(c.isdigit() for c in pw)):
            return pw


def main() -> None:
    parser = argparse.ArgumentParser(description="Provision a SOLANGE login account.")
    parser.add_argument("--email", required=True, help="User's email (their login)")
    parser.add_argument("--name", required=True, help="Full name (shown in the UI)")
    parser.add_argument("--password", default=None,
                        help="Password (omit to auto-generate a strong one)")
    parser.add_argument("--role", default="researcher", choices=["researcher", "admin"],
                        help="Platform role (default: researcher)")
    parser.add_argument("--institution", default="CGU", help="Institution (default: CGU)")
    args = parser.parse_args()

    url = os.environ.get("SUPABASE_URL", "")
    key = os.environ.get("SUPABASE_SERVICE_KEY", "")
    if not url or not key:
        sys.exit("Set SUPABASE_URL and SUPABASE_SERVICE_KEY (service_role) env vars first.")

    sb = create_client(url, key)
    email = args.email.strip().lower()
    password = args.password or generate_password()

    # ── 1. Create the auth user (email auto-confirmed so they can log in now) ──
    user_id = None
    try:
        res = sb.auth.admin.create_user({
            "email": email,
            "password": password,
            "email_confirm": True,
            "user_metadata": {"full_name": args.name},
        })
        user_id = res.user.id
        print(f"✓ Auth account created for {email}")
    except Exception as e:
        msg = str(e)
        if "already" in msg.lower() or "registered" in msg.lower() or "exists" in msg.lower():
            print(f"• Auth account for {email} already exists — looking it up…")
            # Find the existing user id by paging the user list.
            page = 1
            while user_id is None:
                users = sb.auth.admin.list_users(page=page, per_page=200)
                if not users:
                    break
                for u in users:
                    if (u.email or "").lower() == email:
                        user_id = u.id
                        break
                page += 1
            password = None  # unknown / unchanged for existing accounts
        else:
            sys.exit(f"✗ Failed to create auth account: {e}")

    if not user_id:
        sys.exit(f"✗ Could not resolve user id for {email}")

    # ── 2. Ensure the users_profile row exists (name + role) ──────────────────
    try:
        sb.table("users_profile").upsert({
            "id": user_id,
            "full_name": args.name,
            "role": args.role,
            "institution": args.institution,
        }).execute()
        print(f"✓ Profile row set: {args.name} · role={args.role} · {args.institution}")
    except Exception as e:
        print(f"⚠ Auth user is ready but profile row failed: {e}")
        print(f"  You can add it manually in Supabase SQL:")
        print(f"  INSERT INTO users_profile (id, full_name, role, institution) "
              f"VALUES ('{user_id}', '{args.name}', '{args.role}', '{args.institution}');")

    # ── 3. Report credentials to hand to the person ───────────────────────────
    print("\n" + "=" * 52)
    print("  SOLANGE login ready — give these to the user:")
    print("=" * 52)
    print(f"  URL:      https://solange-platform.bio")
    print(f"  Email:    {email}")
    if password:
        print(f"  Password: {password}")
        print("  (They can change it from the in-app account menu.)")
    else:
        print("  Password: unchanged (account already existed)")
    print("=" * 52)


if __name__ == "__main__":
    main()
