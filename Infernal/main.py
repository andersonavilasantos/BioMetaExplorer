import warnings
import subprocess
import os
import time
from Bio import SeqIO
from shutil import copyfile

# Ignora avisos de depreciação
warnings.filterwarnings(action='ignore', category=FutureWarning)
warnings.filterwarnings('ignore')

def remove_equal_sequences(finput, foutput):
    """Remove sequências iguais de um arquivo fasta e grava o resultado em outro arquivo."""
    dataset = {}
    # Abre o arquivo de saída no modo append
    with open(foutput, 'a') as arq:
        for seq_record in SeqIO.parse(finput, "fasta"):
            name = seq_record.name
            seq = seq_record.seq
            if seq not in dataset.values():
                dataset[name] = seq
                arq.write(f">{name}\n{seq}\n")
            else:
                print(f"Removed Equal Sequences: {name}")

start_time = time.time()

print(os.getcwd())
files = os.listdir('Genomes')

for file in files:

	path = "Genomes-Finished/"
	name = file.split('.')[0]
	directory = path + name + '/'

	print(path)
	print(name)
	print(directory)

	if not os.path.exists(directory):
		os.makedirs(directory)

	if os.path.exists('Genomes/' + file):
		copyfile('Genomes/' + file, directory + file)

	subprocess.run(['cmsearch', '--cut_ga', '--rfam', '--nohmmonly','--cpu', '4', '--tblout', str(directory + name + '.tblout'), 'Rfam.cm', str(directory + file)])

	subprocess.run(['esl-sfetch', '--index', str(directory + file)])

	subprocess.run(['sh', 'run_infernal.sh', str(directory + name + '.tblout'), str(directory + file), str(directory + 'seqs.fasta')])

	remove_equal_sequences(str(directory + 'seqs.fasta'), str(directory + 'prep_seqs.fasta'))

end_time = time.time()

execution_time = end_time - start_time

print(f"Execution time: {execution_time} seconds")
