# ── PhenoRewire · Console output ─────────────────────────────────────────────
# What the user sees while a run is going: numbered steps, one or two lines of
# result per step, warnings, and a closing block saying what to open.  Module
# loggers keep logging as before; on screen only the curated progress lines
# and warnings appear, while run_info/run.log keeps everything with timestamps.
#
# Console text is plain ASCII: a Windows console redirected to a file cannot
# encode much beyond it.
from __future__ import annotations

import logging
import sys
import textwrap
from pathlib import Path

PROGRESS_LOGGER = "phenorewire.progress"
progress = logging.getLogger(PROGRESS_LOGGER)

_PACKAGE_LOGGER = "phenorewire"
_WIDTH = 92


class _ConsoleFilter(logging.Filter):
    def __init__(self, verbose: bool) -> None:
        super().__init__()
        self.verbose = verbose

    def filter(self, record: logging.LogRecord) -> bool:
        if self.verbose:
            return True
        return record.name == PROGRESS_LOGGER or record.levelno >= logging.WARNING


class _ConsoleFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        msg = record.getMessage()
        if record.name == PROGRESS_LOGGER and record.levelno < logging.WARNING:
            return msg
        if record.levelno >= logging.ERROR:
            lead = "  ERROR: "
        elif record.levelno >= logging.WARNING:
            lead = "  ! "
        else:
            lead = "  . "
        return "\n".join(_wrap(msg, lead))


def _wrap(msg: str, lead: str) -> list[str]:
    out: list[str] = []
    indent = " " * len(lead)
    for i, para in enumerate(str(msg).splitlines() or [""]):
        first = lead if i == 0 else indent
        out.extend(
            textwrap.wrap(para, width=_WIDTH, initial_indent=first, subsequent_indent=indent)
            or [first.rstrip()]
        )
    return out


def setup_console(level: str = "INFO") -> logging.Handler:
    """Install the console handler on the package logger (once)."""
    pkg = logging.getLogger(_PACKAGE_LOGGER)
    for h in pkg.handlers:
        if getattr(h, "_phenorewire_console", False):
            return h
    # stdout: in a notebook, stderr is shown as if it were an error.
    handler = logging.StreamHandler(sys.stdout)
    handler._phenorewire_console = True  # type: ignore[attr-defined]
    lvl = getattr(logging, str(level).upper(), logging.INFO)
    handler.setLevel(lvl)
    handler.addFilter(_ConsoleFilter(verbose=lvl <= logging.DEBUG))
    handler.setFormatter(_ConsoleFormatter())
    pkg.addHandler(handler)
    pkg.setLevel(logging.DEBUG)
    return handler


def has_any_handler() -> bool:
    return bool(logging.getLogger(_PACKAGE_LOGGER).handlers or logging.getLogger().handlers)


def add_file_log(path: Path) -> logging.Handler:
    """Full DEBUG log of this run, with timestamps, in run_info/run.log."""
    path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(path, mode="w", encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"))
    pkg = logging.getLogger(_PACKAGE_LOGGER)
    pkg.addHandler(handler)
    if pkg.level == logging.NOTSET or pkg.level > logging.DEBUG:
        pkg.setLevel(logging.DEBUG)
    return handler


def remove_handler(handler: logging.Handler) -> None:
    logging.getLogger(_PACKAGE_LOGGER).removeHandler(handler)
    handler.close()


def fmt_p(p) -> str:
    """A p-value as people write it."""
    try:
        p = float(p)
    except (TypeError, ValueError):
        return "NA"
    if p != p:
        return "NA"
    if p < 0.001:
        return "< 0.001"
    return f"{p:.3f}"


def rel(path: Path) -> str:
    """Path relative to the working directory when it is inside it."""
    path = Path(path)
    try:
        return str(path.resolve().relative_to(Path.cwd().resolve())).replace("\\", "/")
    except ValueError:
        return str(path).replace("\\", "/")


class Steps:
    """Numbered progress steps: '[2/6] Preprocessing'."""

    def __init__(self, total: int) -> None:
        self.total = int(total)
        self.done = 0

    def start(self, title: str) -> None:
        self.done += 1
        progress.info("")
        progress.info("[%d/%d] %s", self.done, self.total, title)

    def skip(self, n: int) -> None:
        self.done += int(n)

    @staticmethod
    def detail(text: str, indent: str = "      ") -> None:
        for line in _wrap(text, indent):
            progress.info("%s", line)
