import streamlit as st
from streamlit_option_menu import option_menu
import utils, setup
import dash_bio
from dash import html
import socket
import dash

def is_port_in_use(port: int) -> bool:
    
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) == 0

def runUI():
    st.set_page_config(page_title = "BioMetaExplorer", page_icon = "imgs/icon.png", initial_sidebar_state = "expanded", layout="wide")
    
    utils.inject_css()

    page = option_menu(None, ["Home", "Jobs", "Browse", "About"], 
    icons=["house", "gear-wide", "search", "info-circle"], 
    menu_icon="cast", default_index=0, orientation="horizontal")

    if "queue" not in st.session_state:
        st.session_state["queue"] = False

    if page == "Home":
        setup.home.runUI()
    elif page == "Jobs":
        setup.jobs.runUI()
    elif page == "Browse":
        setup.browse.runUI()
    
    if "dash" not in st.session_state:
        dash_app = dash.Dash(__name__)

        dash_app.layout = html.Div([])

        st.session_state["dash"] = {"app": dash_app}

        for port in range(8050, 8100):
            if not is_port_in_use(port):
                st.session_state["dash"]["port"] = port
                st.session_state["dash"]["app"].run(port=port, debug=False)
                break

if __name__ == "__main__":
    runUI()