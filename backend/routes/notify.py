import os
import smtplib
import logging
from email.message import EmailMessage
from datetime import datetime, timezone

from fastapi import APIRouter, Request

router = APIRouter()
logger = logging.getLogger(__name__)

NOTIFY_TO   = "doron66@gmail.com"
NOTIFY_FROM = "doron66@gmail.com"


def _send_gmail(subject: str, body: str) -> None:
    app_password = os.environ.get("GMAIL_APP_PASSWORD")
    if not app_password:
        logger.warning("GMAIL_APP_PASSWORD not set - login notification skipped")
        return
    msg = EmailMessage()
    msg["Subject"] = subject
    msg["From"]    = NOTIFY_FROM
    msg["To"]      = NOTIFY_TO
    msg.set_content(body)
    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(NOTIFY_FROM, app_password)
        server.send_message(msg)


@router.post("/login")
async def notify_login(request: Request):
    try:
        body       = await request.json()
        user_email = body.get("email", "unknown")
        ip         = request.headers.get("x-forwarded-for", request.client.host if request.client else "unknown")
        ts         = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

        subject = f"SOLANGE login — {user_email}"
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
