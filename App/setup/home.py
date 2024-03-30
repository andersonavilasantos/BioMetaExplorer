import streamlit as st
import polars as pl
import pandas as pd
from io import StringIO
from Bio import SeqIO
import subprocess
import streamlit.components.v1 as components
import os
from secrets import choice
import string
from queue import Queue
from threading import Thread
from streamlit.runtime.scriptrunner import add_script_run_ctx

def submit_job(fasta_sequences, job_path):
    with open(os.path.join(job_path, "predict.fasta"), "w") as f: 
        for record in SeqIO.parse(fasta_sequences, "fasta"):
            f.write(record.format("fasta"))

    subprocess.run(["python", "../Classification/main.py", "--test", job_path] +
                    ["--path_model", "../Classification/results/enc1_cnn_bilstm_4conv_k1_concat2_bio/model.h5"] +
                    ["--encoding", "1", "--k", "1", "--concat", "1", "--feat_extraction", "1", "--features_exist", "1", "--output", job_path])

def queue_listener():
    while True:
        if not job_queue.empty():
            fasta_sequences, job_path = job_queue.get()
            submit_job(fasta_sequences, job_path)
            
def runUI():
    global job_queue

    job_queue = Queue()

    if "listener" not in st.session_state:
        queue_thread = Thread(target=queue_listener)
        add_script_run_ctx(queue_thread)
        queue_thread.start()
        st.session_state["listener"] = True

    img_cols = st.columns([3, 2, 3])

    with img_cols[1]:
        st.image("imgs/logo.png")

    st.markdown("""
        <div style='text-align: center;'>
            <h5 style="color:gray">FAIR and Beyond: Democratizing the Analysis of Non-Coding RNA with Multi-Omics and Data-Centric AI</h5>
        </div>
    """, unsafe_allow_html=True)

    st.info("BioMetaExplorer is ...")

    st.divider()

    st.markdown("""##### Sequence prediction""", unsafe_allow_html=True, help="Only sequences without ambiguous nucleotides are supported.")

    with st.form("sequences_submit"):
        col1, col2 = st.columns(2)

        with col1:
            fasta_text = st.text_area("Input nucleotide sequences in FASTA format", height=125,
                        placeholder=">Sequence_1\nAGCGCAACTCGGACTGCATG\n>Sequence_2\nAGCGGAGTAACTGCATG")
        with col2:
            fasta_file = st.file_uploader("Or upload your FASTA file")

        submitted = st.form_submit_button("Submit")

    predict_path = os.path.abspath("jobs")

    if submitted:
        if fasta_text:
            job_id = ''.join([choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(10)])
            job_path = os.path.join(predict_path, job_id)
            fasta_sequences = StringIO(fasta_text)

            os.makedirs(job_path)
            job_queue.put((fasta_sequences, job_path))
            st.success(f"Job submitted to the queue. You can consult the results in \"Jobs\" using the following ID: **{job_id}**")

        elif fasta_file:
            job_id = ''.join([choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(10)])
            job_path = os.path.join(predict_path, job_id)
            fasta_sequences = StringIO(fasta_file.getvalue().decode("utf-8"))

            os.makedirs(job_path)
            job_queue.put((fasta_sequences, job_path))
            st.success(f"Job submitted to the queue. You can consult the results in \"Jobs\" using the following ID: **{job_id}**")
        else:
            st.error("No sequences submitted!")