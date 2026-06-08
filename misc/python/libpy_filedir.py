
import os
from os.path import expanduser

# Get file size in bytes
def get_filesize(file):
    return os.path.getsize(file)

# Check file existence
def file_exists(file):
    return os.path.exists(file)

# Make directory
def make_directory(dir):
	#	pathlib.Path(dir).mkdir(parents=True, exist_ok=True)
	os.makedirs(dir, exist_ok=True)
	
# Get home directory
def get_homedir():
	return expanduser("~")
