import boto3
from botocore.client import Config
from botocore.utils import fix_s3_host
import urllib3
import argparse
import json
import logging
import sys

# Disable SSL warnings globally as they fill the stderr log otherwise
# Ok to disable because we usually work with local non-443
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def get_sg_s3_connector(credentials_path, logger):
    """
    Get an S3 connector for the specified storage gateway.

    Args:
        credentials_path (str): Path to the credentials file.

    Returns:
        boto3.resource: S3 resource object.

    Raises:
        KeyError: If required configuration keys are missing.
        Exception: For any other exceptions.
    """
    try:
        logger.debug("Initializing S3 connector.")

        # Read credentials from the JSON file
        credentials = json.load(open(credentials_path, 'r'))

        s3_config = Config(signature_version='s3v4',
                           connect_timeout=int(10),  # Time (s) to establish connection
                           read_timeout=int(30),  # Time in seconds to wait for a response
                           retries={'max_attempts': int(20),  # Number of retries
                                    'mode': 'standard'})  # Standard: Exponential backoff, random jitter, 1-8s delay

        # Set up the connection to the local S3 API server
        s3 = boto3.resource('s3',
                            endpoint_url=credentials['endpoint'],
                            aws_access_key_id=credentials['aws_access_key_id'],
                            aws_secret_access_key=credentials['aws_secret_access_key'],
                            verify=False,  # Because we mostly work with local non-443
                            config=s3_config)
        # Not sure if below is even necessary
        s3.meta.client.meta.events.unregister('before-sign.s3', fix_s3_host)

        logger.debug("S3 connector initialized successfully.")
        return s3

    except KeyError:
        logger.error("Missing required configuration key.")
        raise

    except Exception as e:
        logger.error(f"An unexpected error occurred while initializing S3 connector: {e}")
        raise


def setup_logger(name):
    """
    Set up a logger that logs to stdout.

    Args:
        name (str): Name of the logger.

    Returns:
        logging.Logger: Configured logger.
    """
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def download_file(local_path, remote_path, credentials_path, bucket):
    """
    Download a file from an S3 bucket.

    Args:
        local_path (str): Local path to save the downloaded file.
        remote_path (str): Remote path of the file to download.
        credentials_path (str): Path to the credentials file.
        bucket (str): Name of the S3 bucket.

    Raises:
        Exception: If an error occurs during the download.
    """
    logger = setup_logger("hcp_download")

    try:
        logger.debug(f"Downloading file from {remote_path} to {local_path}.")

        # Get the S3 connector
        s3 = get_sg_s3_connector(credentials_path, logger)
        s3_bucket = s3.Bucket(bucket)

        # Download the file
        s3_bucket.download_file(remote_path, local_path)
        logger.debug("File downloaded successfully.")

    except Exception as e:
        logger.error(f"An unexpected error occurred while downloading file: {e}")
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--local_path", help="Local path to save the downloaded file.", required=True, type=str)
    parser.add_argument("-r", "--remote_path", help="Remote path of the file to download.", required=True, type=str)
    parser.add_argument("-c", "--credentials_path", help="Path to the json configuration file.", required=True, type=str)
    parser.add_argument("-b", "--bucket", help="Name of the bucket to download from.", required=True, type=str)
    args = parser.parse_args()

    try:
        download_file(args.local_path, args.remote_path, args.credentials_path, args.bucket)
    except Exception as e:
        print(f"An error occurred: {e}")
        raise


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass