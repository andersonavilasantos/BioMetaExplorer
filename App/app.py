import streamlit as st
from streamlit_option_menu import option_menu
import utils, setup

def runUI():
    st.set_page_config(page_title = "BioMetaExplorer", page_icon = 'imgs/icon.png', initial_sidebar_state = "expanded", layout="wide")
    
    utils.inject_css()

    page = option_menu(None, ["Home", "Browse", "Jobs", "About"], 
    icons=["house", "search", "gear-wide", "info-circle"], 
    menu_icon="cast", default_index=0, orientation="horizontal")

    if page == "Home":
        setup.home.runUI()
    elif page == "Browse":
        setup.browse.runUI()

if __name__ == "__main__":
    runUI()