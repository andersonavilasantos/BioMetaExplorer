import streamlit as st
import os, shutil
import re
import glob
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import polars as pl
import numpy as np

def summary_stats(job_path):
    seq_stats = {"num_seqs": [], "min_length": [], "max_length": [], 
                "avg_length": [], "std_length": [], "sum_length": [], 
                "Q1": [], "Q2": [], "Q3": [], "N50": [], "gc_content": []}
    
    lengths, gcs = [], []
    for record in SeqIO.parse(os.path.join(job_path, "predict.fasta"), "fasta"):
        lengths.append(len(record.seq))
        gcs.append(gc_fraction(record.seq))
    
    num_seqs = len(lengths)

    seq_stats["num_seqs"].append(num_seqs)
    seq_stats["min_length"].append(min(lengths))
    seq_stats["max_length"].append(max(lengths))
    seq_stats["avg_length"].append(np.mean(lengths))
    seq_stats["std_length"].append(np.std(lengths))
    seq_stats["sum_length"].append(sum(lengths))
    seq_stats["Q1"].append(np.percentile(lengths, 25))
    seq_stats["Q2"].append(np.percentile(lengths, 50))
    seq_stats["Q3"].append(np.percentile(lengths, 75))

    lengths.sort(reverse=True)
    total_length = sum(lengths)
    cumulative_length = 0
    for length in lengths:
        cumulative_length += length
        if cumulative_length >= total_length / 2:
            seq_stats["N50"].append(length)
            break
        
    seq_stats["gc_content"].append((sum(gcs) / num_seqs) * 100)
    
    pl.DataFrame(seq_stats).write_csv(os.path.join(job_path, "stats.csv"))