import streamlit as st
from streamlit_option_menu import option_menu
import utils, modules

def runUI():
    st.set_page_config(page_title = "BioMetaExplorer", page_icon = "imgs/icon.png", initial_sidebar_state = "expanded", layout="wide")
    
    utils.inject_css()

    page = option_menu(None, ["Home", "Jobs", "Browse", "About"], 
    icons=["house", "gear-wide", "search", "info-circle"], 
    menu_icon="cast", default_index=0, orientation="horizontal")

    if "queue" not in st.session_state:
        st.session_state["queue"] = False

    if page == "Home":
        modules.home.runUI()
    elif page == "Jobs":
        modules.jobs.runUI()
    elif page == "Browse":
        modules.browse.runUI()

if __name__ == "__main__":
    runUI()