from dash import Dash, dcc, html, dash_table, Input, Output, callback, no_update, State, Patch, DiskcacheManager, clientside_callback
from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash.exceptions import PreventUpdate

import plotly.express as px
import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots
from plotnine import *
import plotnine.options

from PIL import Image
import dash_bootstrap_components as dbc
import scanpy as sc
import os
import pandas as pd
import dask.dataframe as dd
import numpy as np
# import loompy as lp
import h5py
import json

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt
import re
import seaborn as sns
from concurrent import futures

import diskcache
background_callback_manager = DiskcacheManager(diskcache.Cache("/rad/wuc/dash_data/spatial/cache"))

import dash
from dash import Dash, html, dcc

dbc_css = "/home/wuc/dashapps/css/dbc.min.css"
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc_css],
          use_pages=True)

app.layout = html.Div([
    html.H1('Multi-page app with Dash Pages'),
    html.Div([
        html.Div(
            dcc.Link(f"{page['name']} - {page['path']}", href=page["relative_path"])
        ) for page in dash.page_registry.values()
    ]),
    dash.page_container
])

#run server
if __name__ == "__main__":
    app.run(
    host='10.193.0.208',
    port='8051',
    proxy=None,
    debug=False,
    use_reloader=False,
    jupyter_mode='external'
)