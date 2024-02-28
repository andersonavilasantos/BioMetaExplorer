import streamlit as st
import polars as pl
import pandas as pd
from io import StringIO
from Bio import SeqIO
import subprocess
import os

def runUI():
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
    with st.spinner("Predicting sequences..."):
        with st.form("sequences_submit"):
            col1, col2 = st.columns(2)

            with col1:
                fasta_text = st.text_area("Input nucleotide sequences in FASTA format", height=125,
                            placeholder=">Sequence_1\nAGCGCAACTCGGACTGCATG\n>Sequence_2\nAGCGGAGTAACTGCATG")
            with col2:
                fasta_file = st.file_uploader("Or upload your FASTA file")

            submitted = st.form_submit_button("Submit")

        predict_path = os.path.abspath("predict")

        if submitted:

            if fasta_text:
                stringio = StringIO(fasta_text)
                with open(os.path.join(predict_path, "predict.fasta"), "w") as f: 
                    for record in SeqIO.parse(stringio, "fasta"):
                        f.write(record.format("fasta"))

                subprocess.run(["python", "../Classification/main.py", "--test", predict_path] +
                                ["--path_model", "../Classification/results/enc2_cnn_bilstm_4conv_k1_concat1_bio/model.h5"] +
                                ["--encoding", "3", "--k", "1", "--feat_extraction", "1", "--features_exist", "1", "--output", predict_path])
                
                df_results = pd.read_csv("predict/model_predictions.csv")

                df_results["Probability"] = df_results.apply(lambda x: max(x[["Cis-reg", "coding", "rRNA", "sRNA", "tRNA", "unknown"]])*100, axis=1)

                df_results = df_results[["nameseq", "prediction", "Probability"]]
                
                df_results.columns = ["Name", "Prediction", "Probability"]

                st.dataframe(df_results,
                            column_config = {"Probability": st.column_config.ProgressColumn(
                                help="Prediction probability",
                                format="%.2f%%",
                                min_value=0,
                                max_value=100
                            )},
                            hide_index=True, use_container_width=True)
            elif fasta_file:
                stringio = StringIO(fasta_file.getvalue().decode("utf-8"))
                with open(os.path.join(predict_path, "predict.fasta"), "w") as f: 
                    for record in SeqIO.parse(stringio, "fasta"):
                        f.write(record.format("fasta"))

                subprocess.run(["python", "../Classification/main.py", "--test", predict_path] +
                                ["--path_model", "../Classification/results/enc2_cnn_bilstm_4conv_k1_concat1_bio/model.h5"] +
                                ["--encoding", "3", "--k", "1", "--feat_extraction", "1", "--features_exist", "1", "--output", predict_path])
                
                df_results = pl.read_csv("predict/model_predictions.csv")
                
                st.dataframe(df_results.select(["nameseq", "prediction"]), hide_index=True, use_container_width=True)
            else:
                st.error("No sequences submitted!")



