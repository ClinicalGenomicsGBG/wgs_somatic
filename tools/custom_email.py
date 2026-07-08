import os
import smtplib

from definitions import WRAPPER_CONFIG_PATH
from tools.helpers import read_config
from email.message import EmailMessage
from textwrap import dedent, indent

config = read_config(WRAPPER_CONFIG_PATH)
email_config_path = config['email_config_path']
email_config = read_config(email_config_path)
smtp = email_config["smtp_server"]
#mode = "development" if develop else "default" #TODO add develop, right now it defaults to default
mode = "default"
sender = email_config[mode]["sender"]
success_recipients = ", ".join(email_config[mode]["recipients"])
qc_recipients = ", ".join(email_config[mode]["qc"])


def send_email(subject, body):
    """Send a simple email."""

    msg = EmailMessage()
    msg.set_content(body)

    msg["Subject"] = subject
    msg["From"] = sender
    msg["To"] = success_recipients 
    msg["Cc"] = sender

    # Send the message
    s = smtplib.SMTP(smtp)
    s.send_message(msg)
    s.quit()

def send_email_qc(subject, body):
    """Send a simple email."""

    msg = EmailMessage()
    msg.set_content(body)

    msg["Subject"] = subject
    msg["From"] = sender 
    msg["To"] = qc_recipients 
    msg["Cc"] = sender 

    # Send the message
    s = smtplib.SMTP(smtp)
    s.send_message(msg)
    s.quit()


def start_email(run_name, samples):
    """Send an email about starting wgs-somatic for samples in a run"""

    subject = f"WGS Somatic start mail {run_name}"

    sample_list = "\n".join(samples)
    body = f"""\
Starting wgs_somatic for the following samples in run
{run_name}:

{sample_list}

You will get an email when the results are ready.

Best regards,
CGG Cancer
"""

    send_email(subject, body)


def end_email(run_name, samples):
    """Send an email that wgs-somatic has finished running for samples in a run"""

    subject = f"WGS Somatic end mail {run_name}"

    sample_list = "\n".join(samples)
    body = f"""\
WGS somatic has finished successfully for the following samples in run {run_name}:

{sample_list}

Best regards,
CGG Cancer
"""

    send_email(subject, body)


def error_email(run_name, ok_samples, bad_samples):
    """Send an email about which samples have failed and which samples have succeeded"""

    subject = f"Crashed WGS Somatic {run_name}"

    ok_samples_list = "\n".join(ok_samples)
    bad_samples_list = "\n".join(bad_samples)
    body = f"""\
WGS somatic failed for the following samples in run {run_name}:

{bad_samples_list}

The following samples did finish correctly:

{ok_samples_list}

Errors concerning the above samples will be investigated.

Best regards,
CGG Cancer
"""

    send_email(subject, body)


def error_setup_email(instrument):
    """Send an email when the setup of wgs-somatic fails"""

    subject = f"Crashed WGS somatic setup for {instrument}"

    body = f"""\
The automatic setup of WGS somatic failed for instrument {instrument}.

Errors will be investigated.

Best regards,
CGG Cancer
"""

    send_email(subject, body)


def error_admin_qc_email(run_name):
    """Send an email when the generating the qd admin summary report fails"""

    subject = f"WGS somatic - admin QC failed {run_name}"

    body = f"""\
Generating the WGS Admin QC report failed for run {run_name}.

Please create the report manually.
"""

    send_email_qc(subject, body)
