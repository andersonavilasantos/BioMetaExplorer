import os
import subprocess
import time
import warnings
from Bio import SeqIO
from shutil import copyfile
from multiprocessing import Pool

# Ignore deprecation warnings
warnings.filterwarnings(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore')

def remove_equal_sequences(finput, foutput):
    """Remove identical sequences from a fasta file and write the result into another file."""
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

def process_file(file):
    print("Processo~" + file)
    path = "Genomes-Finished"
    name = file.split('.')[0]
    directory = os.path.join(path, name)

    if not os.path.exists(directory):
        os.makedirs(directory)

    source_file = os.path.join('Genomes', file)
    if os.path.exists(source_file):
        dest_file = os.path.join(directory, file)
        copyfile(source_file, dest_file)

    subprocess.run(['cmsearch', '--cut_ga', '--rfam', '--nohmmonly', '--tblout', os.path.join(directory, f"{name}.tblout"), 'Rfam.cm', dest_file])
    subprocess.run(['esl-sfetch', '--index', dest_file])
    subprocess.run(['sh', 'run_infernal.sh', os.path.join(directory, f"{name}.tblout"), dest_file, os.path.join(directory, 'seqs.fasta')])

    remove_equal_sequences(os.path.join(directory, 'seqs.fasta'), os.path.join(directory, 'prep_seqs.fasta'))

files = os.listdir('Genomes')

# Define the number of processes equal to the number of available CPU cores
num_processes = 6

start_time = time.time()

with Pool(processes=num_processes) as pool:
    pool.map(process_file, files)

end_time = time.time()

execution_time = end_time - start_time

print(f"Execution time: {execution_time} seconds")
