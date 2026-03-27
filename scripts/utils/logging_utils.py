#!/usr/bin/env python3
"""
logging_utils.py — Consistent timestamped logging for all Python scripts.
Matches the format of logging_utils.sh so log files read coherently.

Usage:
    from utils.logging_utils import setup_logging, log_info, log_warn, log_error
    from utils.logging_utils import log_step, log_substep, log_separator
    from utils.logging_utils import log_timer_start, log_timer_end, log_table
"""

import sys
import time
import logging
from datetime import datetime
from io import StringIO

# ---------------------------------------------------------------------------
# Module-level state
# ---------------------------------------------------------------------------
_timers: dict[str, float] = {}
_logger: logging.Logger | None = None


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def setup_logging(log_file: str | None = None, level: int = logging.INFO):
    """Configure logging to both stderr and optional file."""
    global _logger
    _logger = logging.getLogger("gwas_pipeline")
    _logger.setLevel(level)
    _logger.handlers.clear()

    fmt = logging.Formatter("%(message)s")

    # Always log to stderr
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    _logger.addHandler(sh)

    # Optionally also log to file
    if log_file:
        fh = logging.FileHandler(log_file, mode="a")
        fh.setFormatter(fmt)
        _logger.addHandler(fh)

    return _logger


def _get_logger() -> logging.Logger:
    global _logger
    if _logger is None:
        setup_logging()
    return _logger


# ---------------------------------------------------------------------------
# Logging functions matching bash equivalents
# ---------------------------------------------------------------------------

def log_info(msg: str):
    _get_logger().info(f"{_timestamp()} [INFO]  {msg}")


def log_warn(msg: str):
    _get_logger().warning(f"{_timestamp()} [WARN]  {msg}")


def log_error(msg: str):
    _get_logger().error(f"{_timestamp()} [ERROR] {msg}")


def log_step(msg: str):
    _get_logger().info(f"\n{_timestamp()} ========== {msg} ==========")


def log_substep(msg: str):
    _get_logger().info(f"{_timestamp()} --- {msg} ---")


def log_separator():
    _get_logger().info(f"{_timestamp()} {'=' * 60}")


def log_timer_start(label: str):
    _timers[label] = time.time()
    log_info(f"Started: {label}")


def log_timer_end(label: str):
    start = _timers.pop(label, time.time())
    elapsed = int(time.time() - start)
    mins, secs = divmod(elapsed, 60)
    log_info(f"Completed: {label} ({mins}m {secs}s)")


def log_table(headers: list[str], rows: list[list], min_width: int = 10):
    """Print a formatted ASCII table to the log."""
    # Calculate column widths
    col_widths = [max(min_width, len(str(h))) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            if i < len(col_widths):
                col_widths[i] = max(col_widths[i], len(str(cell)))

    def _fmt_row(cells):
        parts = []
        for i, cell in enumerate(cells):
            w = col_widths[i] if i < len(col_widths) else min_width
            parts.append(str(cell).ljust(w))
        return "  ".join(parts)

    logger = _get_logger()
    logger.info(f"{_timestamp()} [INFO]  {_fmt_row(headers)}")
    logger.info(f"{_timestamp()} [INFO]  {'  '.join('-' * w for w in col_widths)}")
    for row in rows:
        logger.info(f"{_timestamp()} [INFO]  {_fmt_row(row)}")
