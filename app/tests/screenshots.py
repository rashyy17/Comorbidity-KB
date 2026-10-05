"""Open every page of the running app in headless Chrome, fail on any error element, save screenshots.

Usage (app must be running):  .venv/bin/python app/tests/screenshots.py [base_url]
Uses Playwright with the system Chrome (channel="chrome"); falls back to Playwright's bundled Chromium.
"""
import sys
import time
from pathlib import Path

from playwright.sync_api import sync_playwright

BASE = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:8501"
OUT = Path(__file__).resolve().parents[2] / "results" / "screenshots"
PAGES = [("01_home", ""), ("02_hub_explorer", "hubs"), ("03_gene_profile", "gene"),
         ("04_disease_overlap", "overlap"), ("05_pathways", "pathways"), ("06_network", "network"),
         ("07_expression", "expression"), ("08_validation", "validation"), ("09_query_console", "console")]
SLOW = {"network": 9, "gene": 7}  # pyvis graphs need time for physics to settle


def wait_idle(page, timeout=120):
    """Wait until Streamlit has finished running the script (no 'running' status, no spinners)."""
    t0 = time.time()
    page.wait_for_selector('[data-testid="stMain"]', timeout=timeout * 1000)
    time.sleep(1.5)
    while time.time() - t0 < timeout:
        running = page.locator('[data-testid="stStatusWidget"] :text("Running")').count()
        spinners = page.locator('[data-testid="stSpinner"]').count()
        if not running and not spinners:
            return
        time.sleep(0.5)
    raise TimeoutError("page did not finish rendering")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    problems = []
    with sync_playwright() as p:
        try:
            browser = p.chromium.launch(channel="chrome", headless=True)
        except Exception:
            browser = p.chromium.launch(headless=True)
        ctx = browser.new_context(viewport={"width": 1480, "height": 2400}, device_scale_factor=1.5)
        page = ctx.new_page()
        console_errors = []
        page.on("pageerror", lambda e: console_errors.append(str(e)))
        for name, path in PAGES:
            page.goto(f"{BASE}/{path}")
            wait_idle(page)
            time.sleep(SLOW.get(path, 3))
            exc = page.locator('[data-testid="stException"]').count()
            err = page.locator('[data-testid="stAlertContentError"]').all_inner_texts()
            if exc or err:
                problems.append((path, exc, err))
            # The app scrolls inside stMain: grow the viewport to the content height for a full capture.
            h = page.evaluate("""() => {
                const c = document.querySelector('[data-testid="stMainBlockContainer"]');
                const m = document.querySelector('[data-testid="stMain"]');
                return c ? c.getBoundingClientRect().bottom + (m ? m.scrollTop : 0) : document.body.scrollHeight;
            }""")
            page.set_viewport_size({"width": 1480, "height": min(max(int(h) + 40, 900), 6000)})
            time.sleep(2.5 if path in SLOW else 1.0)
            page.screenshot(path=str(OUT / f"{name}.png"), full_page=True)
            page.set_viewport_size({"width": 1480, "height": 2400})
            print(f"{'ERROR' if exc or err else 'ok':5} {path or 'home':11} -> {name}.png  (height {h}px)")
        browser.close()
    if console_errors:
        print("browser JS errors:", console_errors[:5])
    if problems:
        sys.exit(f"pages with errors: {problems}")
    print("all pages loaded without errors")


if __name__ == "__main__":
    main()
