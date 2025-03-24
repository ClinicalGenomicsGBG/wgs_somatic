import re
import os
import sys


def edit_config(config_template, ploidy, tumor_pileup, normal_pileup, threads, output_config):
    
    outdir = os.path.dirname(output_config)
    os.makedirs(outdir, exist_ok=True)

    with open(config_template, 'r') as file:
        config_data = file.read()
    if ploidy == "free":
        config_data = re.sub(r'^ploidy =.*', 'ploidy = 1,2,3,4,8', config_data, flags=re.MULTILINE)
    else:
        config_data = re.sub(r'^ploidy =.*', f'ploidy = {ploidy}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'^outputDir =.*', f'outputDir = {outdir}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'(\[sample\]\n\nmateFile =).*', f'\\1 {tumor_pileup}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'(\[control\]\n\nmateFile =).*', f'\\1 {normal_pileup}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'^maxThreads =.*', f'maxThreads = {threads}', config_data, flags=re.MULTILINE)
    with open(output_config, 'w') as file:
        file.write(config_data)
    
    return output_config


if __name__ == "__main__":
    config_template = sys.argv[1]
    ploidy = sys.argv[2]
    tumor_pileup = sys.argv[3]
    normal_pileup = sys.argv[4]
    threads = sys.argv[5]
    output_config = sys.argv[6]

    edit_config(config_template, ploidy, tumor_pileup, normal_pileup, threads, output_config)