import os
import json


def find_dropbox_folder():
    # Dropbox's info.json file locations for different operating systems
    paths = {
        "windows": os.path.expandvars(r"%LOCALAPPDATA%\Dropbox\info.json"),
        "mac": os.path.expanduser("~/Library/Application Support/Dropbox/info.json"),
        "linux": os.path.expanduser("~/.dropbox/info.json")
    }

    # Determine the current operating system
    if os.name == 'nt':  # For Windows
        path = paths["windows"]
    elif os.name == 'posix':
        if os.uname().sysname == 'Darwin':  # For macOS
            path = paths["mac"]
        else:  # For Linux
            path = paths["linux"]
    else:
        raise Exception("Unsupported operating system")

    # Check if the info.json file exists
    if os.path.exists(path):
        with open(path, 'r') as file:
            data = json.load(file)
            # The Dropbox folder location can be found in the 'personal' section
            if 'personal' in data and 'path' in data['personal']:
                dropbox_path = data['personal']['path']
                return dropbox_path
            else:
                raise Exception("Could not find the Dropbox folder path in the info.json file")
    else:
        raise Exception("Dropbox info.json file not found")
