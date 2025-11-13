import csv
import logging
import sys

logger = logging.getLogger(__name__)


class SomalierParser:
    """
    Class for parsing a somalier run with a tumor-normal pair or a single tumor or normal sample.

    Attributes:
        relatedness (float): The (first) relatedness value from the pairs file.
        sampleid (str): The sample ID from the normal (if available) or otherwise tumor sample.
        sex (male/female/unknown): The parsed sex prediction from the normal (if available) or otherwise tumor sample.
        y_depth_mean (float): The Y_depth_mean from the normal (if available) or otherwise tumor sample.
        match (bool): Whether the relatedness meets the match cutoff.
    """

    def __init__(
        self,
        pairs_file,
        samples_file,
        tumorstring="tumor",
        normalstring="normal",
        match_cutoff=0.95,
    ):
        self.relatedness = self.parse_somalier_pairs(pairs_file)
        self.sampleid, self.sex, self.y_depth_mean = self.parse_somalier_samples(
            samples_file, tumorstring, normalstring
        )
        if self.relatedness is not None and self.relatedness >= match_cutoff:
            self.match = True
        else:
            self.match = False

    def parse_somalier_pairs(self, pairs_file):
        """
        Parses a somalier pairs.tsv file using csv and returns the relatedness value.
        Assumes data in first row after the header.
        """
        with open(pairs_file, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            row = next(reader, None)
            if row is None:
                logger.info("No relatedness data found in pairs file.")
                return None
            return float(row["relatedness"])

    def parse_somalier_samples(
        self, samples_file, tumorstring="tumor", normalstring="normal"
    ):
        """
        Parses a somalier samples.tsv file using csv.
        Returns the sampleid, sex and Y_depth_mean for the normal if available or otherwise the tumor.
        """
        sex_map = {"1": "male", "2": "female"}
        normal = tumor = None
        with open(samples_file, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                # In the somalier_extract rule, the normal and tumor string are prepended
                if row["sample_id"].startswith(normalstring):
                    normal = (
                        row["sample_id"],
                        sex_map.get(row["sex"], "unknown"),
                        float(row["Y_depth_mean"]),
                    )
                if row["sample_id"].startswith(tumorstring):
                    tumor = (
                        row["sample_id"],
                        sex_map.get(row["sex"], "unknown"),
                        float(row["Y_depth_mean"]),
                    )
        if normal:
            return normal
        elif tumor:
            logger.warning("Normal sample not found; using tumor sample data instead.")
            return tumor
        else:
            logger.error("Neither normal nor tumor data found in samples file.")
            return None, None, None


if __name__ == "__main__":
    # Expect input files as arguments
    pairs_file = sys.argv[1]
    samples_file = sys.argv[2]

    somalier = SomalierParser(pairs_file, samples_file)
    print(somalier.sampleid, somalier.relatedness, somalier.sex)
