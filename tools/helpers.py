import yaml
import json
import logging

def read_config(configpath):
    with open(configpath, 'r') as configfile:
        if configpath.endswith('.json'):
            config_data = json.load(configfile)
        elif configpath.endswith('.yaml') or configpath.endswith('.yml'):
            config_data = yaml.load(configfile, Loader=yaml.FullLoader)
        else:
            raise ValueError("Unsupported file format. Please provide a .json or .yaml file.")
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
