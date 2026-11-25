"""Test original-review design-classification preservation and numerical correction."""
from benchmarks.profile_ahif_review_268 import authors_R_bin,numerically_corrected_bin

def test_source_code_bin_for_18_is_not_equivalent_to_numeric_bin():
    assert authors_R_bin("18")=="No temporal treatment"
    assert numerically_corrected_bin("18")=="11+"

def test_preserve_author_original_review_bins():
    assert authors_R_bin("NA (no temporal treatment)")=="No temporal treatment"
    assert authors_R_bin("1")==authors_R_bin("2")=="1-2"
    assert authors_R_bin("3")==authors_R_bin("10")=="3-10"
    assert authors_R_bin(">10")==">10"

def test_numeric_bin_adjusts_only_missing_18_under_known_review_domain():
    assert numerically_corrected_bin("NA (no temporal treatment)")=="No temporal treatment"
    assert numerically_corrected_bin("1")==numerically_corrected_bin("2")=="1-2"
    assert numerically_corrected_bin("3")==numerically_corrected_bin("10")=="3-10"
    assert numerically_corrected_bin(">10")=="11+"
    assert numerically_corrected_bin("18")=="11+"
