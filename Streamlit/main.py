import streamlit as st
import polars as pl
from streamlit_option_menu import option_menu
import utils
import time

def fetch_data():
    df_pred = pl.read_csv("../Classification/data/predict/results/model_predictions.csv")

    return df_pred

if __name__ == "__main__":
    st.set_page_config(page_title = "BioMetaExplorer", page_icon = 'imgs/icon.png', initial_sidebar_state = "expanded", layout="wide")
    
    utils.inject_css()

    page = option_menu(None, ["Home", "Browse", "About"], 
    icons=["house", "search", "info-circle"], 
    menu_icon="cast", default_index=0, orientation="horizontal")

    if page == "Home":
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

        st.markdown("""##### Sequence prediction</span>
                """, unsafe_allow_html=True, help="Only sequences without ambiguous nucleotides are supported.")

        col1, col2 = st.columns(2)

        with col1:
            st.text_area("Input nucleotide sequences in FASTA format", placeholder=">Sequence_1\nAGCGCAACTCGGACTGCATG")

        with col2:
            st.file_uploader("Or upload your FASTA file")

        st.divider()
    elif page == "Browse":
        # start = time.time()

        st.dataframe(fetch_data())

        # end = time.time()

        # print(end - start)