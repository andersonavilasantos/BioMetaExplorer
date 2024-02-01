import streamlit as st
from streamlit_option_menu import option_menu
import utils

if __name__ == "__main__":
    st.set_page_config(page_title = "BioMetaExplorer", page_icon = 'imgs/icon.png', initial_sidebar_state = "expanded", layout="wide")
    
    utils.inject_css()

    for _ in range(10): 
        st.sidebar.markdown("")
    st.sidebar.divider()

    option = option_menu(None, ["Home", "Browse", "About"], 
    icons=["house", "search", "info-circle"], 
    menu_icon="cast", default_index=0, orientation="horizontal")

    st.markdown('<p style="font-family:Axia Stencil; color:#73b2c1; font-size: 40px;">BioMetaExplorer</p>', unsafe_allow_html=True)

    st.markdown("""##### <span style="color:gray">FAIR and Beyond: Democratizing the Analysis of Non-Coding RNA with Multi-Omics and Data-Centric AI</span>
                """, unsafe_allow_html=True)

    st.info("BioMetaExplorer is ...")

    st.divider()

    st.markdown("""##### Sequence classification</span>
            """, unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        st.text_area("Input nucleotide sequences in FASTA format")

    with col2:
        st.file_uploader("Or upload your FASTA file")

