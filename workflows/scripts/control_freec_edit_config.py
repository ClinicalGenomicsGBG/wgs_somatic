import re
import os
import sys
from sex import calc_sex


def edit_config(config_template, ploidy, tumor_pileup, normal_pileup, chrLenFile, chrFiles, mappability, threads, output_config, wgscovfile, ycovfile):
    # Predict sex using calc_sex
    pred_sex = calc_sex(wgscovfile, ycovfile)

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
    if normal_pileup.lower() == "none":
        config_data = re.sub(r'(\[control\]\n\n)mateFile =.*', r'\1#mateFile =', config_data, flags=re.MULTILINE)
    else:
        config_data = re.sub(r'(\[control\]\n\nmateFile =).*', f'\\1 {normal_pileup}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'^chrLenFile =.*', f'chrLenFile = {chrLenFile}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'^chrFiles =.*', f'chrFiles = {chrFiles}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'^gemMappabilityFile =.*', f'gemMappabilityFile = {mappability}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'^maxThreads =.*', f'maxThreads = {threads}', config_data, flags=re.MULTILINE)
    config_data = re.sub(r'^sex =.*', f'sex = {"XY" if pred_sex == "male" else "XX"}', config_data, flags=re.MULTILINE)
    with open(output_config, 'w') as file:
        file.write(config_data)

    return output_config


if __name__ == "__main__":
    config_template = sys.argv[1]
    ploidy = sys.argv[2]
    tumor_pileup = sys.argv[3]
    normal_pileup = sys.argv[4]
    chrLenFile = sys.argv[5]
    chrFiles = sys.argv[6]
    mappability = sys.argv[7]
    threads = sys.argv[8]
    output_config = sys.argv[9]
    wgscovfile = sys.argv[10]
    ycovfile = sys.argv[11]

    edit_config(config_template, ploidy, tumor_pileup, normal_pileup, chrLenFile, chrFiles, mappability, threads, output_config, wgscovfile, ycovfile)