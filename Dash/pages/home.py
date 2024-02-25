import dash
from dash import html

dash.register_page(__name__, title="BioMetaExplorer", path='/')

layout = html.Div([
    html.Div([
        html.Img(src="assets/logo.png", style={"width": "20%", "display": "block", "margin": "auto"})
    ]),
    
    html.Div([
        html.H1("Your Title"),
        html.P("Your text goes here.")
    ])
])
