import boto3
from botocore.client import Config
from botocore.utils import fix_s3_host
import urllib3
import argparse
import os
from helpers import setup_logger, read_config

# Disable SSL warnings globally as they fill the stderr log otherwise
# Ok to disable because we usually work with local non-443
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)


def get_sg_s3_connector(config, logger):
    """
    Get an S3 connector for the specified storage gateway.

    Args:
        config (ConfigParser): Configuration object containing S3 details.

    Returns:
        boto3.resource: S3 resource object.

    Raises:
        KeyError: If required configuration keys are missing.
        Exception: For any other exceptions.
    """
    try:
        logger.debug("Initializing S3 connector.")

        # Read credentials from the JSON file
        credentials_path = config['hcp']['credentials_file']
        credentials = read_config(credentials_path)

        s3_config = Config(signature_version='s3v4',
                           connect_timeout=int(config['hcp']['connect_timeout']),  # Time (s) to establish connection
                           read_timeout=int(config['hcp']['read_timeout']),  # Time in seconds to wait for a response
                           retries={'max_attempts': int(config['hcp']['max_attempts']),  # Number of retries
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


def download_file(local_path, remote_path, config, logger):
    try:
        logger.debug(f"Downloading file from {remote_path} to {local_path}.")

        # Get the S3 connector
        connector = get_sg_s3_connector(config, logger)
        bucket = connector.Bucket(config['hcp']['bucket'])

        # Download the file
        bucket.download_file(remote_path, local_path)
    
    except Exception as e:
        logger.error(f"An unexpected error occurred while downloading file: {e}")
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--local_path", help="Local path to save the downloaded file.", required=True, type=str)
    parser.add_argument("-r", "--remote_path", help="Remote path of the file to download.", required=True, type=str)
    parser.add_argument("-c", "--config_path", help="Path to the configuration file.", required=True, type=str)
    parser.add_argument("--log_path", help="Path to the log file.", type=str, required=False)
    args = parser.parse_args()


    # Read the configuration file
    config = read_config(args.config_path)

    # Set up the logger
    log_path = args.log_path if args.log_path else config['wrapper_log_path']
    logger = setup_logger("hcp_download", os.path.join(log_path, "hcp_download.log"))

    try:
        download_file(args.local_path, args.remote_path, config, logger)
        logger.info(f"File downloaded successfully to {args.local_path}")

    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        raise


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass