import streamlit as st
import polars as pl
from streamlit_option_menu import option_menu
import plotly.express as px
from math import ceil
import utils
import time

@st.cache_data
def fetch_data():
    df_pred = pl.read_csv("../Classification/data/predict/predictions_complete.tsv", separator="\t")

    return df_pred.to_pandas()

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
            st.text_area("Input nucleotide sequences in FASTA format", 
                        placeholder=">Sequence_1\nAGCGCAACTCGGACTGCATG\n>Sequence_2\nAGCGGAGTAACTGCATG")

        with col2:
            st.file_uploader("Or upload your FASTA file")

        st.divider()
    elif page == "Browse":
        data = fetch_data()

        data_size = len(data)

        table = st.container()

        page_size = 100

        menu = st.columns([6, 1])

        total_pages = ceil(data_size/page_size)

        with menu[1]:
            page_number = st.number_input(
                label="Page",
                label_visibility="collapsed",
                min_value=1,
                max_value=total_pages,
                step=1,
            )

        with menu[0]:
            st.markdown(f"Page **{page_number}** of **{total_pages}** ")

        current_start = (page_number - 1) * page_size

        with table:
            with st.spinner("Loading..."):
                page = data.iloc[current_start:current_start + page_size, :]
                
                col1, col2 = st.columns([2, 1])

                with col1:
                    st.dataframe(page, use_container_width=True)

                with col2:
                    fig = px.sunburst(
                        page,
                        path=["GTDB-tk_domain", "GTDB-tk_phylum", "GTDB-tk_class", "GTDB-tk_order", "GTDB-tk_family", "GTDB-tk_genus", "GTDB-tk_species", "prediction"]
                        # parents="GTDB-tk_domain",
                        # names="GTDB-tk_phylum",
                        
                        # values='nameseq',
                    )
                    
                    st.plotly_chart(fig, use_container_width=True, help="t")