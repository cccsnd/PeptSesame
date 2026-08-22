"""PlantPTM web client — CSRF-authenticated batch submission + polling.

Server protocol (reverse-engineered, see skill references):
    GET  /predict/           -> csrf token in form HTML
    POST /predict/           urlencoded {csrf, fasta_text, ptm_type[], threshold_levels[]}
    GET  /task_result/?task_id=  -> {"completed": bool, "full_data_json": [...]}

Pitfalls encoded here:
- FASTA IDs with `|`-separated metadata make the server return
  completed:True with EMPTY full_data_json, silently — clean IDs only.
- curl `-F field=@file` uploads multipart, Django CharField rejects it;
  must POST via requests `data=` (urlencoded).
"""
from __future__ import annotations

import re
import time
from pathlib import Path

import requests

BASE = "https://ai4bio.online/PlantPTM"
PTM_TYPES = ["Ngly", "Sacy", "Khib", "Kcr", "Ksucc", "Kmal", "Kac", "Kub", "pho"]
THRESHOLDS = ["Extremely high", "High", "Medium", "Low", "Non-PTM"]
UA = {"User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"}


def read_fasta(path: Path) -> list[tuple[str, str]]:
    """Parse FASTA into [(id, seq)] — id is the first whitespace token only."""
    recs: list[tuple[str, str]] = []
    cur_id, cur_seq = None, []
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            if cur_id:
                recs.append((cur_id, "".join(cur_seq)))
            cur_id, cur_seq = line[1:].split()[0], []
        else:
            cur_seq.append(line.strip().upper())
    if cur_id:
        recs.append((cur_id, "".join(cur_seq)))
    return recs


class PlantPTMClient:
    """Thin client for the PlantPTM web service (one request at a time)."""

    def __init__(self, session: requests.Session | None = None):
        self.session = session or requests.Session()
        self.session.headers.update(UA)

    def get_csrf(self) -> str:
        r = self.session.get(f"{BASE}/predict/", timeout=60)
        m = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', r.text)
        if not m:
            raise RuntimeError("csrf token not found")
        return m.group(1)

    def submit_batch(self, fasta: str) -> str:
        """Submit one batch, return task_id."""
        csrf = self.get_csrf()
        data = {
            "csrfmiddlewaretoken": csrf,
            "fasta_text": fasta,
            "ptm_type": PTM_TYPES,
            "threshold_levels": THRESHOLDS,
        }
        r = self.session.post(f"{BASE}/predict/", data=data, timeout=180)
        m = re.search(r'taskId\s*=\s*["\']([^"\']+)["\']', r.text)
        if not m:
            raise RuntimeError(f"no taskId in response (HTTP {r.status_code})")
        return m.group(1)

    def wait_result(self, task_id: str, max_wait: int = 1800) -> dict:
        """Poll task_result until completed; returns result JSON dict."""
        t0 = time.time()
        while time.time() - t0 < max_wait:
            r = self.session.get(f"{BASE}/task_result/?task_id={task_id}", timeout=60)
            d = r.json()
            if d.get("completed"):
                return d
            time.sleep(10)
        raise TimeoutError(f"task {task_id} timeout after {max_wait}s")
