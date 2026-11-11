import pytest

from workflows.scripts import qc_eval


# Run this test once for each set of QC metrics and expected result.
@pytest.mark.parametrize(
    "stype, paired, mean, pct_10x, pct_30x, expected_status",
    [
        pytest.param("tumor", False, 59.9, 0.99, 0.99, "stop", id="tumor-only-below-60"),
        pytest.param("tumor", False, 60, 0.99, 0.99, "pass", id="tumor-only-at-60"),
        pytest.param("tumor", False, 90, 0.99, 0.95, "pass", id="tumor-only-ignores-horizontal"),
        pytest.param("tumor", True, 84, 0.99, 0.95, "stop", id="paired-both-fail"),
        pytest.param("tumor", True, 85, 0.99, 0.95, "warn", id="paired-mean-passes"),
        pytest.param("tumor", True, 84, 0.99, 0.96, "warn", id="paired-horizontal-passes"),
        pytest.param("tumor", True, 85, 0.99, 0.9596, "warn", id="paired-below-96-before-rounding"),
        pytest.param("tumor", True, 85, 0.99, 0.96, "pass", id="paired-both-pass"),
        pytest.param("normal", False, 24, 0.95, 0.99, "stop", id="normal-both-fail"),
        pytest.param("normal", False, 25, 0.95, 0.99, "warn", id="normal-mean-passes"),
        pytest.param("normal", False, 24, 0.96, 0.99, "warn", id="normal-horizontal-passes"),
        pytest.param("normal", False, 25, 0.96, 0.99, "pass", id="normal-both-pass"),
        pytest.param("normal", True, 25, 0.96, 0.99, "pass", id="paired-normal-uses-normal-thresholds"),
    ],
)
def test_qc_status_pass_file_and_email(
    tmp_path, monkeypatch, stype, paired, mean, pct_10x, pct_30x, expected_status
):
    """Coverage metrics determine the status, QC_PASS file, and email subject."""
    # Write the input files needed by evaluate_qc.
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text(
        "# WGS metrics\nMEAN_COVERAGE\tPCT_10X\tPCT_30X\n"
        f"{mean}\t{pct_10x}\t{pct_30x}\n"
    )
    sex = tmp_path / "sex.txt"
    sex.write_text("female\n")
    passfile = tmp_path / "sample.QC_PASS"
    emails = []
    monkeypatch.setattr(qc_eval, "send_email_qc", lambda subject, body: emails.append(subject))

    status, _ = qc_eval.evaluate_coverage(coverage, stype, paired)
    assert status == expected_status

    qc_eval.evaluate_qc(stype, "sample", coverage, sex, passfile, "female", True, paired)
    assert passfile.exists() == (expected_status != "stop")
    expected_emails = {
        "pass": [],
        "warn": ["WGS-somatic QC warning"],
        "stop": ["WGS-somatic has stopped for a sample"],
    }
    assert emails == expected_emails[expected_status]
