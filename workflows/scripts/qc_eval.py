"""
Evaluate WGS QC based on Sentieon coverage metrics and Somalier sex prediction.

Modes:
    Default
            If --coverage is provided, evaluate_coverage is run.
            If --somalier is provided, evaluate_sex is run.

    --evaluate_qc
        evaluate_sex and evaluate_coverage for a sample. Creates a QC_PASS file on success,
        and optionally sends a QC email on fail. 
        Exits with a non-zero status if sample type is not 'tumor' or 'normal'.

Input files:
    --coverage
        Tab-separated output from Sentieon WgsMetricsAlgo.

    --somalier
        Somalier calculated_sex.txt containing the predicted sex.

Outputs:
    - QC_PASS file (if QC passes)
    - Optionally send warning email
"""

import argparse
import sys
import os
from tools.helpers import read_config
from textwrap import dedent
from tools.custom_email import send_email_qc 
from definitions import ROOT_DIR, LAUNCHER_CONFIG_PATH

launcher_config = read_config(LAUNCHER_CONFIG_PATH)
filterconfig = read_config(os.path.join(ROOT_DIR, "configs", launcher_config["filterconf"]))

def evaluate_qc(stype, sname, coverage_file, somalier_file, qc_pass_file, expected_sex=False, email=False, pipestop=False):

    coverage_pass, coverage_message = evaluate_coverage(coverage_file, stype)
    
    if expected_sex:
        somalier_pass, somalier_sex = evaluate_sex(somalier_file, expected_sex)
    else:
        somalier_pass = False
        somalier_sex = "unknown"
    
    message = ''
    if not coverage_pass:
        message += f"""
                    Coverage threshold not reached
                    {coverage_message}
                    ----------------------------------------------------------------------------------------
                    """
        print(dedent(message))
    if not somalier_pass:
        if somalier_sex == "unknown":
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
        print(dedent(message))
    if message:
        if not coverage_pass:
            message = dedent(f"""
                        WGS-somatic pipeline has stopped due to QC error for sample:
                        {sname}
                        ========================================================================================
                        {message}
                        """)
            if email:
                send_email_qc("WGS-somatic has stopped for a sample", message)
                #TODO add mail to geneticists
        elif not somalier_pass:
            message = dedent(f"""
            WGS-somatic is running with a warning for sample:
            {sname}
            ========================================================================================
            {message}
            """)
            if email:
                send_email("QC warning", message)
        else:
            raise ValueError("evaluate_qc crashed due to illogical logic") #Should not happen 

    if coverage_pass:
        print(f"{sname} passed QC")
        with open(qc_pass_file, "w"):
            pass
    else:
        print(f"{sname} failed QC")
    if pipestop:
        print (f"Pipeline will stop")
    else:
        print(f"Pipeline is running")        
        with open(qc_pass_file, "w"):
            pass


def evaluate_coverage(coverage_file, stype):

    with open (coverage_file) as f:
        lines = f.readlines()

    keys = lines[1].strip().split("\t")
    values = [float(x) for x in lines[2].strip().split("\t")]

    wgs_stats = dict(zip(keys, values))
    
    pct_horizontal = filterconfig["qc_pass"]["pct_horizontal"]
    if stype == "tumor":
        cov_threshold = filterconfig["qc_pass"]["tumor_thresholds"]["coverage"] 
        hor_threshold = filterconfig["qc_pass"]["tumor_thresholds"]["horizontal"]
    elif stype == "normal":
        cov_threshold = filterconfig["qc_pass"]["normal_thresholds"]["coverage"]
        hor_threshold = filterconfig["qc_pass"]["normal_thresholds"]["horizontal"]
    else:
        raise ValueError(dedent(f"""\
                Unsupported sample type: {stype}
                Mean coverage: {wgs_stats['MEAN_COVERAGE']} 
                Fraction bases >10X coverage: {wgs_stats['PCT_10X']}
                Fraction bases >30X coverage: {wgs_stats['PCT_30X']}
                """)
                        )

    mean_coverage = wgs_stats['MEAN_COVERAGE']

    try:
        horizontal_coverage = round(wgs_stats[f'PCT_{hor_threshold}X'] * 100, 1)
    except ValueError as e:
        logger.error(dedent(f"""\
        Invalid horizontal coverage value: {e}
        See 'wgs_stats' for possible values
        """))

    message = f"""
                    Mean coverage: {mean_coverage} 
                    Threshold: {cov_threshold}

                    Bases >{hor_threshold}X coverage: {horizontal_coverage}% 
                    Threshold: {pct_horizontal}% of bases >{hor_threshold}X coverage
                    """
    
    if mean_coverage <= cov_threshold or horizontal_coverage <= pct_horizontal:
        return False, message
    return True, message


def evaluate_sex(somalier_file, expected_sex):
    with open(somalier_file) as f:
        sex = f.readline().rstrip()
    
    if sex == expected_sex:
        return True, sex 

    return False, sex 


def main():
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawTextHelpFormatter
            )
    parser.add_argument("-c", "--coverage", help=dedent('''\
                                                        Coverage file. Output from Sentieon WgsMetricsAlgo.
                                                        Expected format:
                                                          - Line 1: Comment beginning with '#'
                                                          - Line 2: Tab-separated column names
                                                          - Line 3: Tab-separated values

                                                        Column names needed: MEAN_COVERAGE, PCT_10X, PCT_30X
                                                        Minimum content of the file, with example values:
                                                        #
                                                        MEAN_COVERAGE    PCT_10X    PCT_30X
                                                        86               98.9       97.7
                                                        ''')
                        )
    parser.add_argument("-s", "--somalier", help="Path to one-line output file from somalier_parse_sex (male/female/unknown)"
    parser.add_argument("-t", "--stype", choices=["tumor", "normal"], default="tumor", help = "Sample type")
    parser.add_argument("-g", "--sex", choices=["male", "female", "unknown"], default = "unknown", help = "Sample sex")
    parser.add_argument("-n", "--name", help = "Sample name")
    parser.add_argument("-e", "--evaluate_qc", action="store_true", help = "Run pipeline behavior: evaluate_qc()")
    parser.add_argument("--email", action="store_true", help = "Send warning email if thresholds are not met")
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    if args.evaluate_qc:
        required = {
        "--coverage": args.coverage,
        "--somalier": args.somalier,
        "--name": args.name
        }
        missing = [arg for arg, value in required.items() if value is None]

        if missing:
            parser.error("--evaluate_qc requires: " + ", ".join(missing))

        evaluate_qc(args.stype, args.name, args.coverage, args.somalier, f"{args.name}.QC_PASS", args.sex, args.email)
       
    else:
        if args.coverage:
            coverage_pass, coverage_metrics = evaluate_coverage(args.coverage, args.stype)
            print(coverage_metrics)

        if args.somalier:
            somalier_pass, somalier_sex = evaluate_sex(args.somalier, args.sex)
            print(f"Expected {args.sex}, Somalier estimate {somalier_sex}. Pass = {somalier_pass}")

if __name__ == "__main__":
    main()
