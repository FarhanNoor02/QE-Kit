import os, json

CONFIG_FILE = os.path.expanduser("~/.qekit_config")

def save_pseudo_dir(path):
    with open(CONFIG_FILE, 'w') as f:
        json.dump({"pseudo_dir": os.path.abspath(path)}, f)

def get_pseudo_dir():
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            return json.load(f).get("pseudo_dir")
    return None
