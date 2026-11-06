import os
import smtplib

from email.message import EmailMessage
from tools.helpers import read_config
from definitions import WRAPPER_CONFIG_PATH

new_line = "\n"

def get_email_settings():
    """Load and return the email settings used for pipeline notifications."""
    config = read_config(WRAPPER_CONFIG_PATH)
    email_config_path = config['email_config_path']
    email_config = read_config(email_config_path)

    smtp = email_config['smtp_server']
    sender = email_config['sender']
    clinic_mail = email_config["clinicians"]
    lab_mail = email_config["lab"]
    bioinfo_mail = email_config["bioinfo"]
    cc = sender

    return smtp, sender, clinic_mail, lab_mail, bioinfo_mail


def send_email(subject, body, clinic=False, lab=False):
    """Send a simple email."""
    smtp, sender, clinic_mail, lab_mail, bioinfo_mail = get_email_settings()

    msg = EmailMessage()
    msg.set_content(body)

    recipients = []

    if not clinic and not lab:
        recipients.extend(bioinfo_mail)
    else:
        if lab:
            recipients.extend(lab_mail)
        if clinic:
            recipients.extend(clinic_mail)
        msg["Cc"] = ", ".join(bioinfo_mail)

    recipients = ", ".join(recipients)

    msg["Subject"] = subject
    msg["From"] = sender
    msg["To"] = recipients

    # Send the message
    with smtplib.SMTP(smtp) as s:
        s.send_message(msg)


def start_email(run_name, samples, clinic=True, lab=True):
    """Send an email about starting wgs-somatic for samples in a run"""

    subject = f"WGS Somatic start mail {run_name}"

    body = f"""Starting wgs_somatic for the following samples in run {run_name}:\n
{new_line.join(samples)}\n
You will get an email when the results are ready.\n
Best regards,
CGG Cancer
"""

    send_email(subject, body, clinic, lab)


def end_email(run_name, ok_samples=[], bad_samples=[], qc_stats=False, clinic=True, lab=True):
    """Send an email that wgs-somatic has finished running for samples in a run"""

    if bad_samples:
        subject = f"Crashed WGS Somatic {run_name}"

        body = f"""WGS somatic failed for the following samples in run {run_name}:\n
{new_line.join(bad_samples)}\n
The following samples did finish correctly:\n
{new_line.join(ok_samples)}\n
Errors concerning the above samples will be investigated.\n
"""

    else:
        subject = f"WGS Somatic end mail {run_name}"

        body = f"""WGS somatic has finished successfully for the following samples in run {run_name}:\n
{new_line.join(ok_samples)}\n
"""
    if not qc_stats:
        body += f"""Generating the WGS Admin QC report failed.\n
Please create the report manually.\n
"""
    
    body +="""
Best regards,
CGG Cancer
"""
    send_email(subject, body, clinic, lab)


def error_setup_email(instrument, clinic=False, lab=True):
    """Send an email when the setup of wgs-somatic fails"""

    subject = f"Crashed WGS somatic setup for {instrument}"

    body = f"""The automatic setup of WGS somatic failed for instrument {instrument}.\n
Errors will be investigated.\n
Best regards,
CGG Cancer
    """

    send_email(subject, body, clinic, lab)
