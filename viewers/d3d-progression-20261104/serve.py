"""Loopback-only static viewer server with deterministic module MIME types."""
from functools import partial
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path

class ViewerHandler(SimpleHTTPRequestHandler):
    extensions_map={**SimpleHTTPRequestHandler.extensions_map,
                    '.mjs':'text/javascript','.js':'text/javascript','.json':'application/json',
                    '.html':'text/html; charset=utf-8','.css':'text/css','.svg':'image/svg+xml'}
    def end_headers(self):
        self.send_header('Cache-Control','no-store')
        self.send_header('X-Content-Type-Options','nosniff')
        super().end_headers()

if __name__=='__main__':
    root=Path(__file__).resolve().parent/'dist'
    ThreadingHTTPServer(('127.0.0.1',8780),partial(ViewerHandler,directory=str(root))).serve_forever()
