#!/usr/bin/env python3
"""
Startup script for the Gene Panel Analysis Server
"""

import os
import sys
import subprocess
from pathlib import Path

def check_requirements():
    """Check if required files exist."""
    required_files = [
        "outputs/deseq2_full_data_no_filter_fc1.5x_p0.05_pvalue.csv",
        "samples.csv",
        "frontend/index.html"
    ]
    
    missing_files = []
    for file_path in required_files:
        if not Path(file_path).exists():
            missing_files.append(file_path)
    
    if missing_files:
        print("❌ Missing required files:")
        for file_path in missing_files:
            print(f"   - {file_path}")
        return False
    
    print("✅ All required files found")
    return True

def install_dependencies():
    """Install Python dependencies."""
    try:
        print("📦 Installing Python dependencies...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"])
        print("✅ Dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed to install dependencies: {e}")
        return False

def start_server():
    """Start the FastAPI server with uvicorn."""
    try:
        print("🚀 Starting Gene Panel Analysis Server (FastAPI)...")
        print("📊 Server will be available at: http://localhost:5001")
        print("📋 API docs will be available at: http://localhost:5001/docs")
        print("🌐 Frontend will be available at: http://localhost:8080")
        print("⚠️  Note: You'll need to serve the frontend folder separately")
        print("   Example: python start_frontend.py")
        print("\n" + "="*60)
        
        # Start uvicorn server
        subprocess.call([
            sys.executable, "-m", "uvicorn", 
            "panel_processor_fastapi:app",
            "--host", "0.0.0.0",
            "--port", "5001",
            "--reload"
        ])
        
    except KeyboardInterrupt:
        print("\n🛑 Server stopped by user")
    except Exception as e:
        print(f"❌ Failed to start server: {e}")

def main():
    """Main startup function."""
    print("🧬 Gene Panel Analysis System Startup")
    print("="*40)
    
    # Check if we're in the right directory
    if not Path("panel_processor_fastapi.py").exists():
        print("❌ Please run this script from the bioinformatics project directory")
        sys.exit(1)
    
    # Check requirements
    if not check_requirements():
        print("\n💡 Please ensure all required data files are in place:")
        print("   - DESeq2 results file in outputs/")
        print("   - samples.csv with patient metadata")
        print("   - Frontend files in frontend/")
        sys.exit(1)
    
    # Install dependencies
    if not install_dependencies():
        sys.exit(1)
    
    # Start server
    start_server()

if __name__ == "__main__":
    main()
