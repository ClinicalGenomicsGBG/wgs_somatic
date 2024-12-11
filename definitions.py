import os

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
WRAPPER_CONFIG_PATH = os.path.join(ROOT_DIR, 'configs', 'wrapper_config.yaml')
LAUNCHER_CONFIG_PATH = os.path.join(ROOT_DIR,'configs','launcher_config.json')
ROOT_LOGGING_PATH = '/clinical/exec/wgs_somatic/logs'