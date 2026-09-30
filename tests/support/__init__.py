"""Test-only helpers shared across the Python test tree; never imported by RelatiPy."""

from pathlib import Path

FIXTURES_DIR = Path(__file__).resolve().parents[1] / "fixtures"
"""Frozen reference data; its paths are recorded in validation provenance."""
