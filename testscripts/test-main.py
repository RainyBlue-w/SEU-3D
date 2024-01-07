from dash import Dash, dcc, dash_table, no_update, State, Patch, DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_daq as daq
import dash_ag_grid as dag
from dash_extensions.enrich import Output, Input, html, callback, DashProxy, LogTransform, DashLogger, DashBlueprint


app = DashProxy(
    __name__, 
    external_stylesheets=[
        dbc.themes.BOOTSTRAP
    ],
    use_pages = True, pages_folder = './pages/',
    external_scripts = [
        {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
    ],
    transforms=[
        LogTransform(), 
    ],
)


app.layout = dbc.Container(
    dbc.Col(
        html.Div()
    )
)
