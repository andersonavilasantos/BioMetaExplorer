import pandas as pd
import subprocess
import argparse

def download_sra(run_access, download_path):
    # Construct the command
    command = f"fastq-dump --split-files --gzip -O {download_path} {run_access}"

    # Execute the command
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Check for any errors
    if result.returncode != 0:
        with open('sra_error.log', 'a') as file:
            file.write(f"Error downloading {run_access}: {result.stderr}\n")
        print(f"Error downloading {run_access}: {result.stderr}")
    else:
        print(f"Download successfully completed: {run_access}")

# Argument parser
parser = argparse.ArgumentParser(description='Download SRA data')
parser.add_argument('--csv', required=True, help='Path to the CSV file')
parser.add_argument('--libsource', required=True, help='Library source to filter, choose between METATRANSCRIPTOMIC or METAGENOMIC')
parser.add_argument('--libstrategy', required=True, help='Library strategy to filter, choose between WGS or AMPLICON')
parser.add_argument('--path', required=True, help='Path to save the downloaded files')
args = parser.parse_args()

# Read the CSV file with pandas
df = pd.read_csv(args.csv)

# Filter dataframe
df = df.loc[(df['Library Source'] == args.libsource) & (df['Library Strategy'] == args.libstrategy)]

# Apply the function to each run_access in the 'run_access' column
for run_access in df['run_access']:
    download_sra(run_access, args.path)
