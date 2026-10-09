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
    msg["Subject"] = subject
    msg["From"] = sender

    recipients = []

    if not clinic and not lab:
        recipients.extend(bioinfo_mail)
    else:
        msg["Cc"] = ", ".join(bioinfo_mail)

        if lab:
            recipients.extend(lab_mail)

        if clinic:
            recipients.extend(clinic_mail)

    msg["To"] = ", ".join(recipients)

    # Send the message
    s = smtplib.SMTP(smtp)
    s.send_message(msg)
    s.quit()


def start_email(run_name, samples):
    """Send an email about starting wgs-somatic for samples in a run"""

    if run_name != "manual":
        subject = f"WGS Somatic start mail {run_name}"

        body = f"""Starting wgs_somatic for the following samples in run {run_name}:\n
{new_line.join(samples)}\n
You will get an email when the results are ready.\n
Best regards,
CGG Cancer
"""
    else:
        subject = "WGS Somatic manual start mail"

        body = f"""Starting wgs_somatic manually for:\n
{new_line.join(samples)}\n
You will get an email when the results are ready.\n
Best regards,
CGG Cancer
"""

    send_email(subject, body, clinic=True, lab=True)


def end_email(run_name, samples):
    """Send an email that wgs-somatic has finished running for samples in a run"""

    if run_name != "manual":
        subject = f"WGS Somatic end mail {run_name}"

        body = f"""WGS somatic has finished successfully for the following samples in run {run_name}:\n
{new_line.join(samples)}\n
Best regards,
CGG Cancer
"""
    else:
        subject = "WGS Somatic manual end mail"

        body = f"""WGS somatic has finished a manual run successfully for the following samples:\n
{new_line.join(samples)}\n
Best regards,
CGG Cancer
"""

    send_email(subject, body, clinic=True, lab=True)


def error_email(run_name, ok_samples=None, bad_samples=None):
    """Send an email about which samples have failed and which samples have succeeded"""

    if run_name != "manual":
        subject = f"Crashed WGS Somatic {run_name}"

        body = f"""WGS somatic failed for the following samples in run {run_name}:\n
{new_line.join(bad_samples)}\n
The following samples did finish correctly:\n
{new_line.join(ok_samples)}\n
Errors concerning the above samples will be investigated.\n
Best regards,
CGG Cancer
"""
    else:
        subject = f"Manual start of WGS Somatic crashed"

        body = f"""WGS somatic failed a manual run for samples:
{new_line.join(bad_samples)}\n

Errors concerning the above samples will be investigated.\n
Best regards,
CGG Cancer
"""

    send_email(subject, body, clinic=True, lab=True)


def error_setup_email(instrument):
    """Send an email when the setup of wgs-somatic fails"""

    subject = f"Crashed WGS somatic setup for {instrument}"

    body = f"""The automatic setup of WGS somatic failed for instrument {instrument}.\n
Errors will be investigated.\n
Best regards,
CGG Cancer
    """

    send_email(subject, body, lab=True)


def error_admin_qc_email(run_name):
    """Send an email when the generating the qd admin summary report fails"""

    subject = f"WGS somatic - admin QC failed {run_name}"

    body = f"""Generating the WGS Admin QC report failed for run {run_name}.\n
Please create the report manually.\n
    """

    send_email(subject, body, lab=True)
