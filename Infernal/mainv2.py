import argparse
import warnings
import subprocess
import os
import time
from Bio import SeqIO
from shutil import copyfile

# Ignore deprecation warnings
warnings.filterwarnings(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore')

def remove_equal_sequences(finput, foutput):
    """Removes identical sequences from a fasta file and writes the result to another file."""
    dataset = {}
    # Open the output file in append mode
    with open(foutput, 'a') as arq:
        for seq_record in SeqIO.parse(finput, "fasta"):
            name = seq_record.name
            seq = seq_record.seq
            if seq not in dataset.values():
                dataset[name] = seq
                arq.write(f">{name}\n{seq}\n")
            else:
                print(f"Removed Equal Sequences: {name}")

# Creating an ArgumentParser object
parser = argparse.ArgumentParser(description="Process the path to the genome folder.")
# Adding argument
parser.add_argument('path', type=str, help="Path to the genome folder.")
# Parse the arguments
args = parser.parse_args()

start_time = time.time()

print(os.getcwd())
files = os.listdir(args.path)  # Use the path argument here

for file in files:

    # Modify path to create new directory in the same location as the original genomes
    path = os.path.join(args.path, 'Genomes-extracted')
    name = file.split('.')[0]
    directory = os.path.join(path, name)

    print(path)
    print(name)
    print(directory)

    if not os.path.exists(directory):
        os.makedirs(directory)

    if os.path.exists(os.path.join(args.path, file)):  # Use the path argument here
        copyfile(os.path.join(args.path, file), os.path.join(directory, file))  # Use the path argument here

    subprocess.run(['cmsearch', '--cut_ga', '--rfam', '--nohmmonly', '--tblout', os.path.join(directory, name + '.tblout'), 'Rfam.cm', os.path.join(directory, file)])

    subprocess.run(['esl-sfetch', '--index', os.path.join(directory, file)])

    subprocess.run(['sh', 'run_infernal.sh', os.path.join(directory, name + '.tblout'), os.path.join(directory, file), os.path.join(directory, 'seqs.fasta')])

    remove_equal_sequences(os.path. join(directory, 'seqs.fasta'), os.path.join(directory, 'prep_seqs.fasta'))

end_time = time.time()

execution_time = end_time - start_time

print(f"Execution time: {execution_time} seconds")


#python3 your_script.py /path/to/directory
