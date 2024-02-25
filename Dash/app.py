from dash import Dash, Input, Output
from dash import html, dcc

import dash_bootstrap_components as dbc

home_layout = html.Div([
    html.Div([
        html.Img(src=r'assets/logo.png', style={'width': '20%', 'display': 'block', 'margin': 'auto'})
    ]),
    
    html.Div([
        html.H1("Your Title"),
        html.P("Your text goes here.")
    ])
])

data_upload_layout = html.Div(
    children=[html.H1(children="This is our upload page")]
)

app = Dash(__name__, title="BioMetaExplorer", external_stylesheets=[dbc.themes.BOOTSTRAP])
# app._favicon = "assets/icon.png"

navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink("Home", href="/")),
        dbc.NavItem(dbc.NavLink("Browse", href="/browse")),
    ],
    brand="BioMetaExplorer",
    color="#73b2c1",
    dark=True,
    className="mb-2",
)

app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        navbar,
        dbc.Container(id="page-content", className="mb-4", fluid=True),
    ]
)


@app.callback(Output("page-content", "children"), Input("url", "pathname"))
def display_page(pathname):
    if pathname == "/":
        return home_layout
    elif pathname == "/browse":
        return data_upload_layout
    else:
        return dbc.Jumbotron(
            [
                html.H1("404: Not found", className="text-danger"),
                html.Hr(),
                html.P(f"The pathname {pathname} was not recognized..."),
            ]
        )


if __name__ == "__main__":
    app.run_server(debug=True)