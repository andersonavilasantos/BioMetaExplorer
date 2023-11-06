import Bio.SeqIO
import pyrodigal
import pyhmmer
from pyhmmer.plan7 import HMMFile
import collections
import subprocess
from tqdm import tqdm
import os
import tempfile

def extract_orfs(genome_files, meta, output):
    for genome in tqdm(genome_files, desc="Extracting ORFs from genomes"):
        with open(output, "a") as file:
            record = Bio.SeqIO.read(genome, "fasta")
            find_orfs(record, meta, file)
        
        remove_duplicates(output)

def find_orfs(record, meta, output_file):
    orf_finder = pyrodigal.GeneFinder(meta=meta)
    if meta:
        genes = orf_finder.find_genes(bytes(record.seq))
    else:
        orf_finder.train(bytes(record.seq))
        genes = orf_finder.find_genes(bytes(record.seq))

    for i, pred in enumerate(genes):
        output_file.write(f">{record.id}_{i+1}\n")
        output_file.write(f"{pred.sequence()}\n")

def remove_duplicates(input_file):
    command = ["seqkit", "rmdup", "-s", input_file]

    completed_process = subprocess.run(command, stdout=subprocess.PIPE, text=True)

    if completed_process.returncode == 0:
        with open(input_file, "w") as output:
            output.write(completed_process.stdout)
    else:
        print("Error: seqkit command failed.")

def translate_sequences(input_file, output_file):
    with open(output_file, "w") as file:
        for record in Bio.SeqIO.parse(input_file, "fasta"):
            file.write(f">{record.id}\n")
            file.write(f"{record.translate().seq}\n")

def identify_pfam_families(input_file, output_file):
    Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])
    
    with pyhmmer.easel.SequenceFile(input_file, digital=True) as seq_file:
        proteins = seq_file.read_block()

    results = []
    with HMMFile("Pfam-A.hmm") as hmm_file:
        for hits in tqdm(pyhmmer.hmmsearch(hmm_file, proteins, bit_cutoffs="trusted"), total=20795, desc="Identifying PFAM families"):
            cog = hits.query_name.decode()
            for hit in hits:
                if hit.included:
                    results.append(Result(hit.name.decode(), cog, hit.score))

    best_results = {}
    keep_query = set()
    for result in results:
        if result.query in best_results:
            previous_bitscore = best_results[result.query].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.query] = result
                keep_query.add(result.query)
            elif result.bitscore == previous_bitscore:
                if best_results[result.query].cog != hit.cog:
                    keep_query.remove(result.query)
        else:
            best_results[result.query] = result
            keep_query.add(result.query)

    filtered_results = [best_results[k] for k in sorted(best_results) if k in keep_query]
    protein_ids = [result.query for result in filtered_results]

    with open(output_file, "w") as file:
        for record in Bio.SeqIO.parse(input_file, "fasta"):
            if record.id in protein_ids:
                file.write(f">{record.id}\n")
                file.write(f"{record.seq}\n")

    return protein_ids

def identify_spurious_proteins(protein_ids, input_file, output_file, original_orfs):
    Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])
    
    with pyhmmer.easel.SequenceFile(input_file, digital=True) as seq_file:
        proteins = seq_file.read_block()

    results = []
    with HMMFile("AntiFam.hmm") as hmm_file:
        for hits in tqdm(pyhmmer.hmmsearch(hmm_file, proteins), total=263, desc="Identifying spurious proteins"):
            cog = hits.query_name.decode()
            for hit in hits:
                if hit.included:
                    results.append(Result(hit.name.decode(), cog, hit.score))

    best_results = {}
    keep_query = set()
    for result in results:
        if result.query in best_results:
            previous_bitscore = best_results[result.query].bitscore
            if result.bitscore > previous_bitscore:
                best_results[result.query] = result
                keep_query.add(result.query)
            elif result.bitscore == previous_bitscore:
                if best_results[result.query].cog != hit.cog:
                    keep_query.remove(result.query)
        else:
            best_results[result.query] = result
            keep_query.add(result.query)

    filtered_results = [best_results[k] for k in sorted(best_results) if k in keep_query]
    spurious_ids = [result.query for result in filtered_results]

    coding_ids = [id for id in protein_ids if id not in spurious_ids]

    with open(output_file, "w") as file:
        for record in Bio.SeqIO.parse(original_orfs, "fasta"):
            if record.id in coding_ids:
                file.write(f">{record.id}\n")
                file.write(f"{record.seq}\n")

if __name__ == "__main__":
    genome_folder = "data"
    meta = False # False if MAG, True if Metagenome
    genome_files = [os.path.join(genome_folder, file) for file in os.listdir(genome_folder) if file.endswith(".fasta")]
    
    # Create a temporary directory to store the files
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_folder = tmp_dir
        orfs_output = os.path.join(tmp_folder, "orfs.fasta")
        proteins_1_output = os.path.join(tmp_folder, "proteins_1.fasta")
        proteins_2_output = os.path.join(tmp_folder, "proteins_2.fasta")
        
        extract_orfs(genome_files, meta, orfs_output)
        translate_sequences(orfs_output, proteins_1_output)
        protein_ids = identify_pfam_families(proteins_1_output, proteins_2_output)
        identify_spurious_proteins(protein_ids, proteins_2_output, "coding.fasta", orfs_output)
        
    print("----- Removed spurious proteins and saved coding sequences -----")