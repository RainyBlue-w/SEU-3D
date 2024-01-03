import dash
from dash import Dash, html, dcc
import dash_bootstrap_components as dbc

dbc_css = "/home/wuc/dashapps/css/dbc.min.css"
app = Dash(
  __name__, 
  external_stylesheets=[
    dbc.themes.BOOTSTRAP,
    "/home/wuc/dashapps/css/dbc.min.css",
  ],
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ],
  use_pages=True
)

header = dbc.NavbarSimple(
    [
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem('spatial', href='/'),
                dbc.DropdownMenuItem("atlas", href='/atlas'),
                dbc.DropdownMenuItem("reik", href='/reik'),
            ],
            nav=True,
            in_navbar=True,
            label="Dataset",
        ),
    ],
    brand="Omics-viewer",
    color="dark",
    dark=True,
    sticky='top',
    style = {"height": "6vh"}
)

app.layout = html.Div([
    header,
    dash.page_container
])



#run server
if __name__ == "__main__":
    app.run(
    host='::',
    port='8050',
    threaded=True,
    proxy=None,
    debug=False,
    use_reloader=False
    # jupyter_mode='external'
)
