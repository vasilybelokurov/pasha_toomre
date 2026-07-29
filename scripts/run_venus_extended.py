#!/usr/bin/env python3
"""Run the reconstructed 13-diagnostic extended Venus document variant."""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from pasha_toomre.cli import main


if __name__ == "__main__":
    raise SystemExit(main(["--planet", "venus", "--layout", "extended", *sys.argv[1:]]))
