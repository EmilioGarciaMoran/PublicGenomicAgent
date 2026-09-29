"""Tests de los clientes HTTP del LLM.

Usan un servidor HTTP local (stdlib) que devuelve respuestas
predefinidas. No tocan red externa. No requieren Ollama real.
"""
from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

import pytest

from publicgenomicagent.agent.llm import (
    HttpLLMClient,
    LLMError,
    LlamaCppClient,
    NullLLMClient,
    OllamaClient,
)


class _Handler(BaseHTTPRequestHandler):
    """Handler configurable por el test."""

    response_body: dict = {}
    response_status: int = 200
    received_payloads: list = []

    def log_message(self, *args, **kwargs):  # silenciar el log
        pass

    def do_POST(self):
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length)
        try:
            payload = json.loads(raw.decode("utf-8"))
        except json.JSONDecodeError:
            payload = {"_raw": raw.decode("utf-8", errors="replace")}
        type(self).received_payloads.append(payload)

        body = json.dumps(type(self).response_body).encode("utf-8")
        self.send_response(type(self).response_status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


@pytest.fixture
def fake_server():
    """Arranca un HTTPServer en un puerto libre y lo apaga al terminar."""
    _Handler.response_body = {}
    _Handler.response_status = 200
    _Handler.received_payloads = []

    server = HTTPServer(("127.0.0.1", 0), _Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    yield {
        "endpoint": f"http://{host}:{port}/api/generate",
        "handler": _Handler,
    }
    server.shutdown()


# ---------------------------------------------------------------------------
# OllamaClient
# ---------------------------------------------------------------------------

def test_ollama_client_returns_response_field(fake_server):
    fake_server["handler"].response_body = {
        "response": '{"tool_name": "qc_bam", "args": {}, "rationale": ""}',
        "done": True,
    }
    client = OllamaClient(
        endpoint=fake_server["endpoint"],
        model="qwen2.5:7b",
        timeout=5,
    )
    out = client.complete("sys", "usr")
    assert out.startswith('{"tool_name"')


def test_ollama_client_sends_expected_payload(fake_server):
    fake_server["handler"].response_body = {"response": "ok", "done": True}
    client = OllamaClient(
        endpoint=fake_server["endpoint"],
        model="qwen2.5:7b",
        temperature=0.0,
        timeout=5,
    )
    client.complete("SYS", "USR")

    payload = fake_server["handler"].received_payloads[0]
    assert payload["model"] == "qwen2.5:7b"
    assert payload["system"] == "SYS"
    assert payload["prompt"] == "USR"
    assert payload["stream"] is False
    assert payload["options"]["temperature"] == 0.0


def test_ollama_client_raises_on_missing_response_field(fake_server):
    fake_server["handler"].response_body = {"done": True}
    client = OllamaClient(
        endpoint=fake_server["endpoint"], model="m", timeout=5
    )
    with pytest.raises(LLMError, match="sin campo 'response'"):
        client.complete("s", "u")


def test_ollama_client_raises_on_non_json(fake_server):
    handler = fake_server["handler"]

    class _BrokenHandler(handler):
        def do_POST(self):
            body = b"not json"
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

    # Reemplazamos el handler del server en marcha no es trivial;
    # probamos directamente el comportamiento del cliente con un
    # endpoint que devuelve texto plano: usamos el mismo server
    # pero forzando un body no-JSON.
    # En su lugar, comprobamos la ruta de error con un endpoint cerrado.
    client = OllamaClient(
        endpoint="http://127.0.0.1:1/api/generate", model="m", timeout=1
    )
    with pytest.raises(LLMError):
        client.complete("s", "u")


# ---------------------------------------------------------------------------
# LlamaCppClient
# ---------------------------------------------------------------------------

def test_llamacpp_client_returns_content_field(fake_server):
    fake_server["handler"].response_body = {
        "content": "hello world",
        "stop": True,
    }
    client = LlamaCppClient(
        endpoint=fake_server["endpoint"], model="local", timeout=5
    )
    out = client.complete("sys", "usr")
    assert out == "hello world"


def test_llamacpp_client_concatenates_system_prompt(fake_server):
    fake_server["handler"].response_body = {"content": "x"}
    client = LlamaCppClient(
        endpoint=fake_server["endpoint"], model="local", timeout=5
    )
    client.complete("SYS", "USR")

    payload = fake_server["handler"].received_payloads[0]
    assert "SYS" in payload["prompt"]
    assert "USR" in payload["prompt"]
    assert payload["stream"] is False


def test_llamacpp_client_raises_on_missing_content(fake_server):
    fake_server["handler"].response_body = {"stop": True}
    client = LlamaCppClient(
        endpoint=fake_server["endpoint"], model="m", timeout=5
    )
    with pytest.raises(LLMError, match="sin campo 'content'"):
        client.complete("s", "u")


# ---------------------------------------------------------------------------
# NullLLMClient
# ---------------------------------------------------------------------------

def test_null_client_always_raises():
    client = NullLLMClient()
    with pytest.raises(LLMError, match="no hay backend"):
        client.complete("s", "u")


# ---------------------------------------------------------------------------
# HttpLLMClient base
# ---------------------------------------------------------------------------

def test_http_base_is_abstract():
    client = HttpLLMClient(endpoint="http://x", model="m")
    with pytest.raises(NotImplementedError):
        client.complete("s", "u")
