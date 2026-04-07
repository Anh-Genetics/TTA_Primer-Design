"""NCBI Entrez API setup and efetch wrapper with retry / rate-limiting.

Public functions:
    setup_entrez   — configure Biopython Entrez (email, api_key, rate limit).
    entrez_efetch  — call Entrez.efetch with rate limiting and exponential
                     back-off retry.
"""

from __future__ import annotations

import logging
import time
from typing import IO

from tta_primer_design.utils.rate_limiter import RateLimiter

logger = logging.getLogger("tta_primer_design.utils.ncbi_api")

# Module-level limiter shared across all callers in the same process.
_limiter: RateLimiter | None = None


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def setup_entrez(email: str, api_key: str = "") -> None:
    """Configure Biopython Entrez and initialise the rate limiter.

    NCBI allows **3 requests/second** without an API key and **10/second**
    with one.  Call this once at application start-up (or whenever credentials
    change).

    Args:
        email:   Required by NCBI; used to identify the caller.
        api_key: Optional NCBI API key.  When provided the rate limit is
                 raised from 3 to 10 requests/second.
    """
    global _limiter
    from Bio import Entrez  # lazy import — keeps startup fast if Bio unavailable

    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
        rate = 10.0
    else:
        # Clear any previously set key so the slower limit is applied.
        if hasattr(Entrez, "api_key"):
            Entrez.api_key = ""
        rate = 3.0

    _limiter = RateLimiter(calls_per_second=rate)
    logger.debug("Entrez configured: email=%s, api_key_set=%s", email, bool(api_key))


def entrez_efetch(
    *,
    db: str,
    id: str,
    rettype: str,
    retmode: str = "text",
    retries: int = 3,
    timeout: int = 120,
) -> IO[str]:
    """Call ``Entrez.efetch`` with rate limiting and exponential-backoff retry.

    The module-level rate limiter is applied before every attempt.  On
    transient errors (network timeouts, HTTP 429, etc.) the call is retried
    up to *retries* times with an increasing sleep between attempts
    (2^0, 2^1, … seconds).

    Args:
        db:      Entrez database name (e.g. ``"nucleotide"``).
        id:      Accession or numeric UID to fetch.
        rettype: Return type (``"fasta"``, ``"gb"``, ``"gbc"``, …).
        retmode: Return mode (default ``"text"``).
        retries: Maximum number of attempts before raising.
        timeout: Per-request timeout in seconds (informational; passed
                 through for future Entrez HTTP timeout support).

    Returns:
        File-like handle returned by ``Bio.Entrez.efetch``.

    Raises:
        RuntimeError: When all *retries* attempts fail.
    """
    from Bio import Entrez

    if _limiter is None:
        # Auto-configure with whatever email is already set (if any).
        setup_entrez(getattr(Entrez, "email", None) or "user@example.com")

    assert _limiter is not None  # satisfied by setup_entrez above

    last_exc: Exception | None = None
    for attempt in range(retries):
        try:
            _limiter.acquire()
            handle = Entrez.efetch(db=db, id=id, rettype=rettype, retmode=retmode)
            return handle
        except Exception as exc:
            last_exc = exc
            wait = 2**attempt
            logger.warning(
                "Entrez.efetch attempt %d/%d failed for id=%r: %s. Retrying in %ds…",
                attempt + 1,
                retries,
                id,
                exc,
                wait,
            )
            time.sleep(wait)

    raise RuntimeError(
        f"Entrez.efetch failed after {retries} attempts for id={id!r}: {last_exc}"
    ) from last_exc
