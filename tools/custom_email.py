import os
import smtplib

from email.message import EmailMessage

new_line = '\n'

def send_email(subject, body):
    """Send a simple email."""

    msg = EmailMessage()
    msg.set_content(body)

    msg['Subject'] = subject
    msg['From'] = "cgg-cancer@gu.se" # TODO Get from config
    msg['To'] = "gms_btb@gu.se" # TODO Get from config and have different recipients for errors and success
    msg['Cc'] = "cgg-cancer@gu.se" # TODO Get from config


    # Send the message
    s = smtplib.SMTP('smtp.gu.se')
    s.send_message(msg)
    s.quit()


def start_email(run_name, samples):
    """Send an email about starting wgs-somatic for samples in a run"""

    subject = f'WGS Somatic start mail {run_name}'

    body = f"""Starting wgs_somatic for the following samples in run {run_name}:\n
{new_line.join(samples)}\n
You will get an email when the results are ready.\n
Best regards,
CGG Cancer
 """

    send_email(subject, body)

def end_email(run_name, samples):
    """Send an email that wgs-somatic has finished running for samples in a run"""

    subject = f'WGS Somatic end mail {run_name}'

    body = f"""WGS somatic has finished successfully for the following samples in run {run_name}:\n
{new_line.join(samples)}\n
Best regards,
CGG Cancer
"""

    send_email(subject, body)


def error_email(run_name, ok_samples, bad_samples):
    """Send an email about which samples have failed and which samples have succeeded"""

    subject = f'WGS Somatic end mail {run_name}'

    body = f"""WGS somatic has finished successfully for the following samples in run {run_name}:\n
{new_line.join(ok_samples)}\n
The following samples have not finished correctly:\n
{new_line.join(bad_samples)}\n
Errors concerning these samples will be investigated.\n
Best regards,
CGG Cancer
"""

    send_email(subject, body)
