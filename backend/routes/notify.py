import os
import smtplib
import logging
import traceback
from email.message import EmailMessage
from datetime import datetime, timezone

from fastapi import APIRouter, Request

router = APIRouter()
logger = logging.getLogger(__name__)

NOTIFY_TO   = "doron66@gmail.com"
NOTIFY_FROM = "doron66@gmail.com"


def _send_gmail(subject: str, body: str) -> None:
    app_password = os.environ.get("GMAIL_APP_PASSWORD", "")
    if not app_password:
        logger.warning("GMAIL_APP_PASSWORD not set - login notification skipped")
        return

    # Diagnose any non-ASCII in the password
    if not app_password.isascii():
        non_ascii = [(i, repr(c)) for i, c in enumerate(app_password) if ord(c) > 127]
        logger.error("GMAIL_APP_PASSWORD contains non-ASCII chars: %s", non_ascii)
        return

    # Sanitize all inputs — strip anything non-ASCII before it hits SMTP
    subject = subject.encode("ascii", errors="replace").decode("ascii")
    body    = body.encode("ascii", errors="replace").decode("ascii")

    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"]    = NOTIFY_FROM
    msg["To"]      = NOTIFY_TO
    msg.set_content(body)

    try:
        with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
            server.login(NOTIFY_FROM, app_password)
            server.send_message(msg)
        logger.info("Login notification sent")
    except Exception:
        logger.error("SMTP send failed:\n%s", traceback.format_exc())


@router.post("/login")
async def notify_login(request: Request):
    try:
        body       = await request.json()
        user_email = body.get("email", "unknown")
        ip         = request.headers.get("x-forwarded-for", request.client.host if request.client else "unknown")
        ts         = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

        subject = f"SOLANGE login - {user_email}"
        text    = (
            f"SOLANGE platform login\n\n"
            f"User:  {user_email}\n"
            f"Time:  {ts}\n"
            f"IP:    {ip}\n"
        )
        _send_gmail(subject, text)
    except Exception as e:
        logger.error("Login notification failed: %s", e)

    return {"ok": True}
