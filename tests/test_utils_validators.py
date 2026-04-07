"""Tests cho utils/validators.py."""

from __future__ import annotations

import pytest

from tta_primer_design.utils.validators import validate_accession, validate_email, validate_sequence

# ---------------------------------------------------------------------------
# validate_accession
# ---------------------------------------------------------------------------


class TestValidateAccession:
    """Test NCBI accession format validator."""

    # --- valid RefSeq ---
    @pytest.mark.parametrize(
        "acc",
        [
            "NM_001101",
            "NM_001101.5",
            "NR_024540",
            "NG_007400",
            "XM_017021954",
            "XR_001234567",
            "NP_000001",
            "NC_000001",
            "NT_000001",
            "NW_000001",
            "NZ_000001000",
            "WP_000000001",
            "AC_000001",
            "AP_000001",
        ],
    )
    def test_valid_refseq(self, acc: str) -> None:
        assert validate_accession(acc) is True

    # --- valid GenBank ---
    @pytest.mark.parametrize(
        "acc",
        [
            "U12345",
            "AY123456",
            "AAB12345",
            "U12345.1",
            "AY123456.2",
        ],
    )
    def test_valid_genbank(self, acc: str) -> None:
        assert validate_accession(acc) is True

    # --- invalid ---
    @pytest.mark.parametrize(
        "acc",
        [
            "",
            "   ",
            "not_an_accession",
            "12345",
            "NM",
            "NM_",
            "NM_abc",
            "TOOLONG_123456789012",
        ],
    )
    def test_invalid(self, acc: str) -> None:
        assert validate_accession(acc) is False

    def test_none_like_empty_string(self) -> None:
        assert validate_accession("") is False

    def test_non_string_returns_false(self) -> None:
        assert validate_accession(None) is False  # type: ignore[arg-type]

    def test_strips_whitespace(self) -> None:
        assert validate_accession("  NM_001101  ") is True


# ---------------------------------------------------------------------------
# validate_sequence
# ---------------------------------------------------------------------------


class TestValidateSequence:
    """Test nucleotide sequence validator."""

    @pytest.mark.parametrize("seq", ["ATCG", "atcg", "ATCGN", "AAUUCG", "NNNN"])
    def test_valid_standard(self, seq: str) -> None:
        assert validate_sequence(seq) is True

    @pytest.mark.parametrize("seq", ["ATCGR", "ATCGY", "ATCGW"])
    def test_iupac_rejected_by_default(self, seq: str) -> None:
        assert validate_sequence(seq) is False

    @pytest.mark.parametrize("seq", ["ATCGR", "ATCGY", "RYSWKMBDHV"])
    def test_iupac_accepted_when_allowed(self, seq: str) -> None:
        assert validate_sequence(seq, allow_iupac=True) is True

    def test_empty_string(self) -> None:
        assert validate_sequence("") is False

    def test_non_string(self) -> None:
        assert validate_sequence(None) is False  # type: ignore[arg-type]

    def test_invalid_characters(self) -> None:
        assert validate_sequence("ATCG123") is False
        assert validate_sequence("ATCG ") is False

    def test_mixed_case_valid(self) -> None:
        assert validate_sequence("AtCgAtCg") is True

    def test_all_n(self) -> None:
        assert validate_sequence("NNNNNN") is True


# ---------------------------------------------------------------------------
# validate_email
# ---------------------------------------------------------------------------


class TestValidateEmail:
    """Test email validator."""

    @pytest.mark.parametrize(
        "email",
        [
            "user@example.com",
            "researcher@ncbi.nlm.nih.gov",
            "a.b+c@domain.org",
            "test123@test.co.uk",
        ],
    )
    def test_valid_emails(self, email: str) -> None:
        assert validate_email(email) is True

    @pytest.mark.parametrize(
        "email",
        [
            "",
            "notanemail",
            "@domain.com",
            "user@",
            "user@domain",
            "user @domain.com",
        ],
    )
    def test_invalid_emails(self, email: str) -> None:
        assert validate_email(email) is False

    def test_non_string(self) -> None:
        assert validate_email(None) is False  # type: ignore[arg-type]

    def test_strips_whitespace(self) -> None:
        assert validate_email("  user@example.com  ") is True
