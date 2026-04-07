"""Tests cho utils/rate_limiter.py."""

from __future__ import annotations

import threading
import time
from unittest.mock import patch

import pytest

from tta_primer_design.utils.rate_limiter import RateLimiter


class TestRateLimiterInit:
    """Test RateLimiter construction."""

    def test_default_rate(self) -> None:
        limiter = RateLimiter()
        assert limiter._min_interval == pytest.approx(1.0 / 3.0)

    def test_custom_rate(self) -> None:
        limiter = RateLimiter(calls_per_second=10.0)
        assert limiter._min_interval == pytest.approx(0.1)

    def test_invalid_rate_raises(self) -> None:
        with pytest.raises(ValueError, match="calls_per_second"):
            RateLimiter(calls_per_second=0)

    def test_negative_rate_raises(self) -> None:
        with pytest.raises(ValueError):
            RateLimiter(calls_per_second=-1.0)


class TestRateLimiterAcquire:
    """Test acquire() behaviour."""

    def test_first_acquire_does_not_sleep(self) -> None:
        limiter = RateLimiter(calls_per_second=3.0)
        with patch("time.sleep") as mock_sleep:
            limiter.acquire()
            mock_sleep.assert_not_called()

    def test_second_acquire_sleeps_when_called_too_fast(self) -> None:
        limiter = RateLimiter(calls_per_second=1.0)  # 1 s interval
        slept: list[float] = []

        def capture_sleep(secs: float) -> None:
            slept.append(secs)
            # Don't actually sleep in the test — just record the value.

        with patch("time.sleep", side_effect=capture_sleep):
            # Force _last_call_time to be "just now" so the next acquire will
            # see that insufficient time has elapsed.
            limiter._last_call_time = time.monotonic()
            limiter.acquire()

        assert len(slept) == 1
        assert slept[0] > 0

    def test_acquire_allows_call_after_interval(self) -> None:
        limiter = RateLimiter(calls_per_second=100.0)
        # Set last call far in the past — acquire should not sleep.
        limiter._last_call_time = time.monotonic() - 10.0
        with patch("time.sleep") as mock_sleep:
            limiter.acquire()
            mock_sleep.assert_not_called()

    def test_last_call_time_updated(self) -> None:
        limiter = RateLimiter(calls_per_second=3.0)
        before = time.monotonic()
        limiter.acquire()
        after = time.monotonic()
        assert before <= limiter._last_call_time <= after


class TestRateLimiterContextManager:
    """Test context-manager protocol."""

    def test_context_manager_enter_returns_self(self) -> None:
        limiter = RateLimiter(calls_per_second=100.0)
        with limiter as ctx:
            assert ctx is limiter

    def test_context_manager_calls_acquire(self) -> None:
        limiter = RateLimiter(calls_per_second=100.0)
        acquired: list[bool] = []
        original_acquire = limiter.acquire

        def track_acquire() -> None:
            acquired.append(True)
            original_acquire()

        limiter.acquire = track_acquire  # type: ignore[method-assign]
        with limiter:
            pass
        assert acquired == [True]


class TestRateLimiterThreadSafety:
    """Basic thread-safety smoke test."""

    def test_concurrent_acquires_do_not_raise(self) -> None:
        limiter = RateLimiter(calls_per_second=1000.0)
        errors: list[Exception] = []

        def worker() -> None:
            try:
                limiter.acquire()
            except Exception as exc:
                errors.append(exc)

        threads = [threading.Thread(target=worker) for _ in range(10)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        assert errors == []
