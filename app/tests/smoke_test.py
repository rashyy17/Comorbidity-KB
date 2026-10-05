"""Headless smoke test: render every page with Streamlit's AppTest and fail on any exception.

Run from the repository root:  .venv/bin/python app/tests/smoke_test.py
"""
import sys
from pathlib import Path

from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parent.parent
PAGES = ["home", "hubs", "gene", "overlap", "pathways", "network", "expression", "validation", "console"]

SCRIPT = """
import sys
sys.path.insert(0, {app!r})
import streamlit as st
from views import {page}
{page}.render()
"""


def main():
    failed = []
    for page in PAGES:
        at = AppTest.from_string(SCRIPT.format(app=str(APP), page=page), default_timeout=180)
        at.run()
        errors = [e.value for e in at.exception] + [e.body for e in at.error]
        status = "OK" if not errors else "FAIL"
        print(f"{status:4} {page:11} widgets={len(at.main.children):3} {errors[:2] if errors else ''}")
        if errors:
            failed.append(page)
    if failed:
        sys.exit(f"failed pages: {failed}")
    print("all pages rendered without errors")


if __name__ == "__main__":
    main()
