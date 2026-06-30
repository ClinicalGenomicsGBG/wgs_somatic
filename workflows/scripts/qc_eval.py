import argparse
import sys
from tools.custom_email import send_email_qc 

def evaluate_qc(stype, sname, coverage_file, somalier_file, qc_pass_file, expected_sex=False):

    coverage_pass, coverage_message = evaluate_coverage(coverage_file, stype)
    
    if expected_sex:
        somalier_pass, somalier_sex = evaluate_sex(somalier_file, expected_sex)
    else:
        somalier_pass = False
        somalier_sex = "Missing"
    
    message = ''
    if not coverage_pass:
        message += f"""
Coverage threshold not reached
{coverage_message}
----------------------------------------------------------------------------------------
"""
    if not somalier_pass:
        if somalier_sex == "Missing":
            message += """
Somalier warning:
No information regarding sex for patient. 
"""
        else:
            message += f"""
Somalier warning:
Patient is: {expected_sex}
Somalier estimate: {somalier_sex}
""" 

    if message:
        message = f"""
WGS-somatic pipeline is running with warnings for sample:
{sname}
========================================================================================
{message}
"""
        send_email_qc("QC warning", message)


    with open(qc_pass_file, "w"):
        pass

def evaluate_coverage(coverage_file, stype):
    with open (coverage_file) as f:
        lines = f.readlines()

    keys = lines[1].strip().split("\t")
    values = [float(x) for x in lines[2].strip().split("\t")]

    wgs_stats = dict(zip(keys, values))
    if stype == "tumor":
        horizontal_coverage = wgs_stats['PCT_30X']
        cov_threshold = 85
        hor_threshold = 30
    elif stype == "normal":
        cov_threshold = 25
        hor_threshold = 10
    else:
        raise ValueError(f"""Unsupported sample type: {stype}
Median coverage: {wgs_stats['MEDIAN_COVERAGE']} 
Fraction bases >10X coverage: {wgs_stats['PCT_10X']}
Fraction bases >30X coverage: {wgs_stats['PCT_30X']}
""")

    median_coverage = wgs_stats['MEDIAN_COVERAGE']
    horizontal_coverage = round(wgs_stats[f'PCT_{hor_threshold}X'] * 100, 1)

    message = f"""Median coverage: {median_coverage} 
Threshold: {cov_threshold}

Bases >{hor_threshold}X coverage: {horizontal_coverage}% 
Threshold: 95% of bases >{hor_threshold}X coverage
"""
    
    if median_coverage <= cov_threshold or horizontal_coverage <= 95:
        return False, message
    return True, message

def evaluate_sex(somalier_file, expected_sex):
    with open(somalier_file) as f:
        sex = f.readline().rstrip()
    
    if sex == expected_sex:
        return True, sex 
    return False, sex 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--coverage")
    parser.add_argument("-s", "--somalier")
    parser.add_argument("-t", "--stype", choices=["tumor", "normal"], default="tumor")
    parser.add_argument("-g", "--sex", choices=["male", "female"])
    parser.add_argument("-n", "--name")
    parser.add_argument("-e", "--evaluate_qc", action="store_true", help = "Run evaluate_qc() to call coverage and sex functions as in pipeline")
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if args.evaluate_qc:
        try:
            evaluate_qc(args.stype, args.name, args.coverage, args.somalier, "/dev/null" , args.sex)
        except:
            parser.parse_args(['--help'])
    else:
        if args.coverage:
            coverage_pass, coverage_metrics = evaluate_coverage(args.coverage, args.stype)
            print(coverage_metrics)

        if args.somalier:
            somalier_pass, somalier_sex = evaluate_sex(args.somalier, args.sex)
            print(f"Expected {args.sex}, Somalier estimate {somalier_sex}. Pass = {somalier_pass}")

if __name__ == "__main__":
    main()
