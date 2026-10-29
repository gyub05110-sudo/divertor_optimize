"""Local-only viewer server with deterministic JavaScript MIME types on Windows."""
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
import argparse
class ViewerHandler(SimpleHTTPRequestHandler):
    extensions_map={**SimpleHTTPRequestHandler.extensions_map,'.mjs':'text/javascript','.js':'text/javascript',
                    '.json':'application/json','.css':'text/css','.html':'text/html'}
    def __init__(self,*args,**kwargs):
        super().__init__(*args,directory=str(Path(__file__).resolve().parent),**kwargs)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8765)
    args=parser.parse_args()
    ThreadingHTTPServer(('127.0.0.1',args.port),ViewerHandler).serve_forever()
