from dash import Dash, dcc, dash_table, no_update, State, Patch, DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import dash_daq as daq
import dash_ag_grid as dag
from dash_extensions.enrich import Output, Input, html, callback, DashProxy, LogTransform, DashLogger, DashBlueprint

import pandas as pd

layout = dbc.Col(
    dbc.Input(id='input'),
    dbc.Button(id = 'button')
)