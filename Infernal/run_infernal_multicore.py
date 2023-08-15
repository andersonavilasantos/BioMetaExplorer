
import argparse
import os
import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor

def run_mainv2(input_path, output_path):
    """
    Execute mainv2.py as a subprocess.
    
    Parameters:
    - input_path (str): Path to the genome folder.
    - output_path (str): Output path for the extracted folder.
    """
    try:
        script_path = 'mainv2.py'
        result = subprocess.run(['python3', script_path, '-p', input_path, '-o', output_path], capture_output=True, text=True)
        print(result.stdout)
        if result.stderr:
            print(f"Errors:\n{result.stderr}")

    except Exception as e:
        print(f"An error occurred: {e}")


def process_folder(folder):
    parent_folder = os.path.dirname(folder.rstrip(os.sep))
    infernal_folder = os.path.join(parent_folder, 'infernal')

    if not os.path.exists(infernal_folder):
        os.makedirs(infernal_folder)
        print(f"'infernal' folder created at {infernal_folder}")

    run_mainv2(folder, infernal_folder)


def list_folders(path):
    """Lists the folders based on the provided glob pattern."""
    try:
        folders = glob.glob(path)
        
        if not folders:
            print(f"No folders found for the path: {path}")
            return

        # Using ThreadPoolExecutor to run processes in parallel
        with ThreadPoolExecutor(max_workers=6) as executor:
            executor.map(process_folder, folders)

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List folders based on the provided glob path pattern.")
    parser.add_argument('-p', '--path', required=True, help='Glob path pattern to specify which folders to list. E.g. "/path/to/folder/*/"')
    args = parser.parse_args()
    list_folders(args.path)
