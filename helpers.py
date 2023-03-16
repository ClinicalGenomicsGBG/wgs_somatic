import yaml
import json
import logging
from definitions import ROOT_DIR

def read_config(configpath):
    with open(configpath, 'r') as configfile:
        config_data = json.load(configfile)
        return config_data

def read_clusterconf():
    with open(f"{ROOT_DIR}/configs/cluster.yaml") as inputfile:
        inputfile_info = yaml.load(inputfile, Loader=yaml.FullLoader)
        return inputfile_info

def read_passconfig():
    with open("/root/password_config.json", 'r') as configfile:
        config_data = json.load(configfile)
        return config_data

def setup_logger(name, log_path=None):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    if log_path:
        file_handle = logging.FileHandler(log_path, 'a')
        file_handle.setLevel(logging.DEBUG)
        file_handle.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(module)s - %(message)s'))
        logger.addHandler(file_handle)

    return logger

