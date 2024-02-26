import streamlit as st
import polars as pl
from math import ceil
import plotly.express as px

@st.cache_data(show_spinner=False)
def fetch_data():
    df_pred = pl.read_csv("../Classification/data/predict/predictions_complete.tsv", separator="\t")

    return df_pred.to_pandas()

def runUI():
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
                st.dataframe(page, height=500, hide_index=True, use_container_width=True)

            with col2:
                fig = px.sunburst(
                    page,
                    path=["GTDB-tk_domain", "GTDB-tk_phylum", "GTDB-tk_class", "GTDB-tk_order", "GTDB-tk_family", "GTDB-tk_genus", "GTDB-tk_species", "prediction"]
                    # parents="GTDB-tk_domain",
                    # names="GTDB-tk_phylum",
                    
                    # values='nameseq',
                )
                
                st.plotly_chart(fig, use_container_width=True, help="t")