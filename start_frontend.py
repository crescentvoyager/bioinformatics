#!/usr/bin/env python3
"""
Simple HTTP server for the frontend
"""

import http.server
import socketserver
import os
import webbrowser
from pathlib import Path

PORT = 8080
FRONTEND_DIR = "frontend"

def main():
    """Start the frontend server."""
    
    # Check if frontend directory exists
    if not Path(FRONTEND_DIR).exists():
        print(f"❌ Frontend directory '{FRONTEND_DIR}' not found")
        return
    
    # Change to frontend directory
    os.chdir(FRONTEND_DIR)
    
    # Create server
    Handler = http.server.SimpleHTTPRequestHandler
    httpd = socketserver.TCPServer(("", PORT), Handler)
    
    print("🌐 Gene Panel Analysis Frontend")
    print("="*35)
    print(f"🚀 Server starting on http://localhost:{PORT}")
    print("📝 Make sure the backend is running on http://localhost:5001")
    print(f"📁 Serving files from: {Path.cwd()}")
    print("\nPress Ctrl+C to stop the server")
    print("="*35)
    
    # Open browser
    try:
        webbrowser.open(f"http://localhost:{PORT}")
    except:
        pass
    
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\n🛑 Frontend server stopped")
        httpd.shutdown()

if __name__ == "__main__":
    main()
