"""fetch() del script de calibración: reintento, reanudación y comprobación."""
from __future__ import annotations

import importlib.util
import os
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

SCRIPT = Path(__file__).parents[2] / "scripts" / "calibrate_relatedness_giab.py"
DATA = os.urandom(300_000)


def _load():
    spec = importlib.util.spec_from_file_location("calib_fetch", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


class _Flaky(BaseHTTPRequestHandler):
    """Primera petición: anuncia todo y corta a mitad. Después: atiende Range."""
    calls: list[str] = []

    def do_GET(self):
        rng = self.headers.get("Range")
        _Flaky.calls.append(rng or "full")
        if rng:
            start = int(rng.split("=")[1].split("-")[0])
            body = DATA[start:]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(DATA) - 1}/{len(DATA)}")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        self.send_response(200)
        self.send_header("Content-Length", str(len(DATA)))
        self.end_headers()
        self.wfile.write(DATA[: len(DATA) // 2])     # corte a mitad
        self.wfile.flush()
        self.connection.shutdown(2)

    def log_message(self, *a):
        pass


@pytest.fixture
def server():
    _Flaky.calls = []
    srv = HTTPServer(("127.0.0.1", 0), _Flaky)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_port}/ALL.chr22.vcf.gz"
    srv.shutdown()


def test_fetch_resumes_after_cut(server, tmp_path):
    calib = _load()
    out = calib.fetch(server, tmp_path / "f.vcf.gz")
    assert out.read_bytes() == DATA
    assert _Flaky.calls[0] == "full" and _Flaky.calls[1].startswith("bytes=")
    assert not (tmp_path / "f.vcf.gz.part").exists()


def test_fetch_integrity_check_failure_aborts(server, tmp_path):
    calib = _load()
    with pytest.raises(SystemExit):
        calib.fetch(server, tmp_path / "g.vcf.gz", retries=2, check=lambda p: False)
    assert not (tmp_path / "g.vcf.gz").exists()


def test_fetch_is_cached(server, tmp_path):
    calib = _load()
    dest = tmp_path / "h.bin"
    dest.write_bytes(b"ya estaba")
    assert calib.fetch(server, dest).read_bytes() == b"ya estaba"
    assert _Flaky.calls == []
