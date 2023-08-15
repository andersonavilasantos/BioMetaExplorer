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

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def list_folders(path):
    """Lists the folders based on the provided glob pattern."""
    try:
        folders = glob.glob(path)
        
        if not folders:
            print(f"No folders found for the path: {path}")
            return

        # Limit to 18 processes
        max_processes = 12

        # Split the folders list into chunks of 18 folders each
        folder_chunks = list(chunks(folders, max_processes))

        for folder_chunk in folder_chunks:
            # Using ThreadPoolExecutor to run processes in parallel
            with ThreadPoolExecutor(max_workers=max_processes) as executor:
                executor.map(process_folder, folder_chunk)

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="List folders based on the provided glob path pattern.")
    parser.add_argument('-p', '--path', required=True, help='Glob path pattern to specify which folders to list. E.g. "/path/to/folder/*/"')
    args = parser.parse_args()
    list_folders(args.path)
