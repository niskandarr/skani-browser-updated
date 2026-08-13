from http.server import HTTPServer, SimpleHTTPRequestHandler

class COIHandler(SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        super().end_headers()

if __name__ == "__main__":
    server = HTTPServer(("0.0.0.0", 8080), COIHandler)
    print("Serving on port 8080")
    server.serve_forever()