#!/usr/bin/env python3
"""Generate LEON's Ed25519 signing keypair (P8 seal non-repudiation upgrade).

Run this ONCE, then:
  1. Set LEON_SIGNING_PRIVATE_KEY on the backend deployment (Render) to the
     printed private key -- an environment secret, NEVER committed to the repo.
  2. Set LEON_PUBLIC_KEY the same way (needed by routes.provenance's reverify
     path to check a stored signature) -- this one is also safe to publish
     (dissertation, docs): it is a public key, not a secret.

Losing the private key does not invalidate any existing SHA-256 integrity
digest (p8_hash) -- those remain valid and checkable. It only means runs
sealed from that point on are signed with a NEW identity; anything signed
under the old key can no longer be verified against the new public key,
so treat this as a one-time setup, not something to regenerate casually.
"""
import base64
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey

if __name__ == "__main__":
    key = Ed25519PrivateKey.generate()
    private_b64 = base64.b64encode(key.private_bytes_raw()).decode()
    public_b64 = base64.b64encode(key.public_key().public_bytes_raw()).decode()
    print("LEON_SIGNING_PRIVATE_KEY (secret -- set on Render, never commit):")
    print(f"  {private_b64}")
    print()
    print("LEON_PUBLIC_KEY (safe to publish -- dissertation, docs, other env vars):")
    print(f"  {public_b64}")
