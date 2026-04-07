"""Input validators for TTA Primer Design.

Functions:
    validate_accession  — check NCBI accession format (RefSeq & classic GenBank).
    validate_sequence   — check nucleotide sequence characters.
    validate_email      — basic syntactic email check required by Entrez.
"""

from __future__ import annotations

import re

# ---------------------------------------------------------------------------
# Compiled patterns
# ---------------------------------------------------------------------------

# RefSeq prefixes: NC, NG, NM, NR, NT, NW, NZ, NP, XM, XP, XR, WP, AC, AP
_REFSEQ_PATTERN = re.compile(
    r"^(?:NC|NG|NM|NR|NT|NW|NZ|NP|XM|XP|XR|WP|AC|AP)_\d{6,9}(?:\.\d+)?$",
    re.IGNORECASE,
)

# Classic GenBank: 1–2 letters + 5–6 digits  (e.g. U12345, AY123456)
#                  3 letters + 5 digits        (e.g. AAB12345)
#                  4–6 letters + 8–10 digits   (WGS contigs / proteins)
_GENBANK_PATTERN = re.compile(
    r"^(?:" r"[A-Z]{1,2}\d{5,6}" r"|[A-Z]{3}\d{5,8}" r"|[A-Z]{4,6}\d{8,10}" r")(?:\.\d+)?$",
    re.IGNORECASE,
)

# Valid nucleotide characters (strict IUPAC)
_VALID_BASES = frozenset("ACGTUNacgtun")
_VALID_IUPAC = frozenset("ACGTURYSWKMBDHVNacgturyswkmbdhvn")

# Minimal email: something@something.something
_EMAIL_PATTERN = re.compile(r"^[^@\s]+@[^@\s]+\.[^@\s]+$")


# ---------------------------------------------------------------------------
# Public validators
# ---------------------------------------------------------------------------


def validate_accession(accession: str) -> bool:
    """Check whether *accession* matches a known NCBI accession format.

    Accepts RefSeq accessions (e.g. ``NM_001101``, ``NM_001101.5``) and
    classic GenBank accessions (e.g. ``U12345``, ``AY123456``).

    Args:
        accession: String to validate.

    Returns:
        ``True`` if the format matches a known pattern, ``False`` otherwise.
    """
    if not accession or not isinstance(accession, str):
        return False
    acc = accession.strip()
    return bool(_REFSEQ_PATTERN.match(acc) or _GENBANK_PATTERN.match(acc))


def validate_sequence(sequence: str, *, allow_iupac: bool = False) -> bool:
    """Check whether *sequence* is a valid nucleotide string.

    Args:
        sequence: Nucleotide sequence to validate.
        allow_iupac: When ``True``, accept full IUPAC ambiguity codes
                     (R, Y, S, W, K, M, B, D, H, V) in addition to the
                     standard bases A, C, G, T, U, N.

    Returns:
        ``True`` if every character is a valid nucleotide symbol.
    """
    if not sequence or not isinstance(sequence, str):
        return False
    valid = _VALID_IUPAC if allow_iupac else _VALID_BASES
    return all(c in valid for c in sequence)


def validate_email(email: str) -> bool:
    """Check whether *email* is syntactically valid for use with NCBI Entrez.

    Only basic format validation is performed (presence of ``@`` and a
    domain with at least one dot).  Full RFC-5322 compliance is **not**
    guaranteed.

    Args:
        email: Email address string to validate.

    Returns:
        ``True`` if the format appears valid, ``False`` otherwise.
    """
    if not email or not isinstance(email, str):
        return False
    return bool(_EMAIL_PATTERN.match(email.strip()))
