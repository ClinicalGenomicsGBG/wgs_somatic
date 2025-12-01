import yaml
import json
import logging
import os
import glob


def read_config(configpath):
    with open(configpath, "r") as configfile:
        if configpath.endswith(".json"):
            config_data = json.load(configfile)
        elif configpath.endswith(".yaml") or configpath.endswith(".yml"):
            config_data = yaml.load(configfile, Loader=yaml.FullLoader)
        else:
            raise ValueError(
                "Unsupported file format. Please provide a .json or .yaml file."
            )
        return config_data


def setup_logger(name, log_path=None):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    if log_path:
        file_handle = logging.FileHandler(log_path, "a")
        file_handle.setLevel(logging.DEBUG)
        file_handle.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(module)s - %(message)s")
        )
        logger.addHandler(file_handle)

    return logger


def collect_versions(version_dir, outpath, extension="*.txt"):
    """
    Collect version information from all files with the given extension in version_dir

    version_dir: Directory containing version files
    outpath: Path to output the collected versions
    extension: File extension to look for (default: .txt)
    """
    dirpath = os.path.dirname(outpath)
    if dirpath:
        os.makedirs(dirpath, exist_ok=True)

    files = sorted(glob.glob(os.path.join(version_dir, extension)))

    try:
        with open(outpath, "w") as out:
            out.write("# Collected Versions\n\n")
            if not files:
                out.write(
                    f"# No version files found in {version_dir} matching {extension}.\n"
                )
                return

            for f in files:
                key = os.path.splitext(os.path.basename(f))[
                    0
                ]  # filename without extensions as heading
                out.write(f"{key}: |\n")

                try:
                    with open(f, "r") as infile:
                        for line in infile:
                            out.write(
                                f"    {line}"
                            )  # Add indentation for YAML block style
                except Exception as e_file:
                    out.write(f"    # Error reading file {f}: {e_file}\n")

                out.write("\n")  # Add a newline between tools
    except Exception as e:
        raise RuntimeError(
            f"Failed to write collected versions to {outpath}: {e}"
        ) from e
