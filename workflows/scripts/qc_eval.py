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
from typing import Tuple
from tools.helpers import read_config
from textwrap import dedent
from tools.custom_email import send_email 
from definitions import ROOT_DIR, LAUNCHER_CONFIG_PATH

launcher_config = read_config(LAUNCHER_CONFIG_PATH)
filterconfig = read_config(os.path.join(ROOT_DIR, "configs", launcher_config["filterconf"]))

def evaluate_qc(stype, sname, coverage_file, somalier_file, qc_pass_file, expected_sex=False, send_email=False, paired=False):

    coverage_status, coverage_message = evaluate_coverage(coverage_file, stype, paired)
    
    if expected_sex:
        somalier_pass, somalier_sex = evaluate_sex(somalier_file, expected_sex)
    else:
        somalier_pass = False
        somalier_sex = "unknown"
    
    message = ''
    if coverage_status != "pass":
        message += f"Coverage threshold not reached\n{coverage_message}{'-' * 88}\n"
    if not somalier_pass:
        if somalier_sex == "unknown":
            message += "Somalier warning:\nNo information regarding sex for patient.\n"
        else:
            message += f"Somalier warning:\nPatient sex is: {expected_sex}\nSomalier estimate: {somalier_sex}\n"
        print(message)
    if message:
        if coverage_status == "stop":
            message = f"WGS-somatic pipeline has stopped due to QC error for sample:\n{sname}\n{'=' * 88}\n{message}"
            subject = "WGS-somatic has stopped for a sample"
        else:
            message = f"WGS-somatic is running with a warning for sample:\n{sname}\n{'=' * 88}\n{message}"
            subject = "WGS-somatic QC warning"
        if send_email:
            send_email_qc(subject, message)

    if coverage_status != "stop":
        print(f"{sname} passed QC")
        with open(qc_pass_file, "w"):
            pass
    else:
        print(f"{sname} failed QC", file=sys.stderr)
        print(message, file=sys.stderr)
        print(f"{sname} failed QC", file=sys.stdout)
        print(message, file=sys.stdout)


def evaluate_coverage(coverage_file, stype, paired=False) -> Tuple[str, str]:
    """Return (status, message), with status "pass", "warn", or "stop".

    With two configured checks, one failure warns and two failures stop.
    With one check, failure stops.
    The message contains measured values and limits.
    """

    with open (coverage_file) as f:
        lines = f.readlines()

    keys = lines[1].strip().split("\t")
    values = [float(x) for x in lines[2].strip().split("\t")]

    wgs_stats = dict(zip(keys, values))
    
    pct_horizontal = filterconfig["qc_pass"]["pct_horizontal"]
    if stype == "tumor":
        threshold_group = "tumor_paired_thresholds" if paired else "tumor_only_thresholds"
        tumor_thresholds = filterconfig["qc_pass"][threshold_group]
        cov_threshold = tumor_thresholds["coverage"]
        hor_threshold = tumor_thresholds["horizontal"]
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

    horizontal_coverage = None
    if hor_threshold is not None:
        horizontal_coverage = wgs_stats[f'PCT_{hor_threshold}X'] * 100

    mean_ok = cov_threshold is None or mean_coverage >= cov_threshold
    horizontal_ok = hor_threshold is None or horizontal_coverage >= pct_horizontal

    message = f"Mean coverage: {mean_coverage}\n"
    if cov_threshold is not None:
        message += f"Threshold: {cov_threshold}X\n"
    if hor_threshold is not None:
        message += f"\nBases >{hor_threshold}X coverage: {horizontal_coverage:.1f}%\n"
        message += f"Threshold: {pct_horizontal}% of bases >{hor_threshold}X coverage\n"

    if cov_threshold is None and hor_threshold is None:
        return "pass", message
    if cov_threshold is None:
        return ("pass" if horizontal_ok else "stop"), message
    if hor_threshold is None:
        return ("pass" if mean_ok else "stop"), message
    if mean_ok and horizontal_ok:
        return "pass", message
    if mean_ok or horizontal_ok:
        return "warn", message
    return "stop", message


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
    parser.add_argument("-s", "--somalier", help="Path to one-line output file from somalier_parse_sex (male/female/unknown)")
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
            _, coverage_metrics = evaluate_coverage(args.coverage, args.stype)
            print(coverage_metrics)

        if args.somalier:
            somalier_pass, somalier_sex = evaluate_sex(args.somalier, args.sex)
            print(f"Expected {args.sex}, Somalier estimate {somalier_sex}. Pass = {somalier_pass}")

if __name__ == "__main__":
    main()
