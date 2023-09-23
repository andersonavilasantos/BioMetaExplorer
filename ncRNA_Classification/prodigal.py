from Bio import SeqIO
import pyrodigal
import os

# Define the input and output folders
input_folder = "Infernal/genomes_train"
output_folder = "test"

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# Iterate through all .fasta files in the input folder
for filename in os.listdir(input_folder):
    if filename.endswith(".fasta"):
        input_filepath = os.path.join(input_folder, filename)
        output_filepath = os.path.join(output_folder, filename)

        # Read the input FASTA file
        record = SeqIO.read(input_filepath, "fasta")

        # Create an ORF finder
        orf_finder = pyrodigal.GeneFinder()  # Use "meta" if it's a metagenome
        orf_finder.train(bytes(record.seq))
        genes = orf_finder.find_genes(bytes(record.seq))

        # Write the genes to the output file
        with open(output_filepath, "w") as dst:
            genes.write_genes(dst, sequence_id=record.id)