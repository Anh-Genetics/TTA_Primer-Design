"""Thread-safe token-bucket rate limiter for NCBI API calls."""

from __future__ import annotations

import threading
import time


class RateLimiter:
    """Simple interval-based rate limiter (thread-safe).

    Ensures that at most ``calls_per_second`` calls are made per second by
    sleeping for the remainder of the minimum inter-call interval when
    ``acquire()`` is invoked more frequently than allowed.

    Args:
        calls_per_second: Maximum allowed calls per second (must be > 0).

    Example::

        limiter = RateLimiter(calls_per_second=3)
        with limiter:
            response = requests.get(url)
    """

    def __init__(self, calls_per_second: float = 3.0) -> None:
        if calls_per_second <= 0:
            raise ValueError(f"calls_per_second must be > 0, got {calls_per_second}")
        self._min_interval: float = 1.0 / calls_per_second
        self._lock = threading.Lock()
        self._last_call_time: float = 0.0

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def acquire(self) -> None:
        """Block the calling thread until the next call is permitted."""
        with self._lock:
            now = time.monotonic()
            wait = self._min_interval - (now - self._last_call_time)
            if wait > 0:
                time.sleep(wait)
            self._last_call_time = time.monotonic()

    # ------------------------------------------------------------------
    # Context-manager support
    # ------------------------------------------------------------------

    def __enter__(self) -> RateLimiter:
        self.acquire()
        return self

    def __exit__(self, *args: object) -> None:
        pass
