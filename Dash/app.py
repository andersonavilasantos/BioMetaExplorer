import dash
from dash import Dash, Input, Output
from dash import html, dcc

import dash_bootstrap_components as dbc

app = Dash(__name__, 
            title="BioMetaExplorer", 
            external_stylesheets=[dbc.themes.BOOTSTRAP],
            use_pages=True)

navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink("Home", href="/", active="exact")),
        dbc.NavItem(dbc.NavLink("Browse", href="/browse", active="exact")),
    ],
    brand="BioMetaExplorer",
    color="#73b2c1",
    dark=True,
    className="mb-2",
)

app.layout = html.Div(
    [
        navbar,
        dash.page_container
    ]
)

if __name__ == "__main__":
    app.run_server(debug=True)