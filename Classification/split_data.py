from Bio import SeqIO
import random
import os

def split_fasta_in_directory(directory, train_ratio=0.8):
    train_dir = os.path.join(directory, "train")
    test_dir = os.path.join(directory, "test")

    if not os.path.exists(train_dir):
        os.mkdir(train_dir)
    if not os.path.exists(test_dir):
        os.mkdir(test_dir)

    fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]

    for fasta_file in fasta_files:
        sequences = list(SeqIO.parse(os.path.join(directory, fasta_file), "fasta"))

        random.shuffle(sequences)

        num_train = int(len(sequences) * train_ratio)
        train_sequences = sequences[:num_train]
        test_sequences = sequences[num_train:]

        SeqIO.write(train_sequences, os.path.join(train_dir, fasta_file), "fasta")
        SeqIO.write(test_sequences, os.path.join(test_dir, fasta_file), "fasta")

# Uso
directory_path = "/mnt/UFZ-Data/anderson/NAR-Data"
split_fasta_in_directory(directory_path)
