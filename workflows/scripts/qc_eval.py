import argparse

def evaluate_qc(stype, coverage_file, somalier_file, qc_pass_file):

    coverage_pass, message = evaluate_coverage(coverage_file, stype)
    
    expected_sex = "male" #placeholder for now
    somalier_pass, somalier_sex = evaluate_sex(somalier_file, expected_sex)

    if not coverage_pass and somalier_pass:

        errors = []

        if not cov_ok:
            errors.append(f"Coverage failed: {message}")

        if not som_ok:
            errors.append(f"Somalier failed: {somalier_sex}")

        if errors:
            raise ValueError("\n".join(errors))

        with open(qc_pass_file, "w"):
            pass

def evaluate_coverage(coverage_file, stype):
    with open (coverage_file) as f:
        lines = f.readlines()

    keys = lines[1].strip().split("\t")
    values = [float(x) for x in lines[2].strip().split("\t")]

    wgs_stats = dict(zip(keys, values))
    
    if stype == "tumor":
        message = f"Median coverage: {wgs_stats["MEDIAN_COVERAGE"]}, Fraction bases with >30X coverage: {wgs_stats["PCT_30X"]}"
        if wgs_stats["MEDIAN_COVERAGE"] <= 85 or wgs_stats["PCT_30X"] <= 0.95:
            return False, message
        else:
            return True, message

    elif stype == "normal":
        message = f"Median coverage: {wgs_stats["MEDIAN_COVERAGE"]}, Fraction bases with >10X coverage: {wgs_stats["PCT_10X"]}"
        if wgs_stats["MEDIAN_COVERAGE"] <= 25 or wgs_stats["PCT_10X"] <= 0.95:
            return False, message 
        else:
            return True, message

    else:
        raise ValueError("Unupported sample type")

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
    args = parser.parse_args()

    if args.coverage:
        coverage_pass, coverage_metrics = evaluate_coverage(args.coverage, args.stype)
        print(coverage_metrics)

    if args.somalier:
        somalier_pass, somalier_sex = evaluate_sex(args.somalier, args.sex)
        print(f"Expected {args.sex}, Somalier estimate {somalier_sex}. Pass = {somalier_pass}")

if __name__ == "__main__":
    main()
