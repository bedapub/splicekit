import os
import splicekit
import http.server
import socket
import socketserver
import RangeHTTPServer

def process():
    setup()
    server()

def server():
    port = 8005
    hostname=socket.gethostname()
    ip_addr=socket.gethostbyname(hostname)
    Handler = RangeHTTPServer.RangeRequestHandler
    socketserver.TCPServer.allow_reuse_address = True
    with socketserver.TCPServer(("", port), Handler) as httpd:
        print(f"[splicekit] http://{ip_addr}:{port}")
        httpd.serve_forever()

def setup():
    core_path = os.path.dirname(splicekit.__file__)
    jbrowse2_folder = "jbrowse2"
    if not os.path.exists(jbrowse2_folder):
        os.makedirs(jbrowse2_folder)
    if not os.path.exists("jbrowse2/version.txt"):
        os.system("wget https://github.com/GMOD/jbrowse-components/releases/download/v2.4.2/jbrowse-web-v2.4.2.zip -O jbrowse2/jbrowse-web-v2.4.2.zip")
        os.system("unzip -qq -d jbrowse2 -o jbrowse2/jbrowse-web-v2.4.2.zip")
        os.system("rm jbrowse2/jbrowse-web-v2.4.2.zip")
