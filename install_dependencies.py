import os
import platform
import subprocess

def main():
    cmd = "pip install --upgrade -r requirements.txt --user"
    
    if platform.system() == "Windows":
        # Execute in PowerShell
        subprocess.run(["powershell", "-Command", cmd], check=True)
    elif platform.system() == "Linux":
        # Execute in Shell
        subprocess.run(["/bin/sh", "-c", cmd], check=True)
    else:
        print("Unsupported OS")

if __name__ == "__main__":
    main()
