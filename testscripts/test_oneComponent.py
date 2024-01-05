# In[] env
from math import isnan
import math
from dash import Dash, dcc, html, dash_table, Input, Output, callback, no_update, State, Patch, DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash.exceptions import PreventUpdate
import dash_daq as daq
import dash_ag_grid as dag

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

import squidpy as sq

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt
import re
import seaborn as sns
from concurrent import futures
from typing import List, Dict, Tuple
import diskcache
background_callback_manager = DiskcacheManager(diskcache.Cache("/rad/wuc/dash_data/spatial/cache"))


# In[] data

exp_data = {
  'E7.5': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.5_HC0.5_min400.h5ad"),
  'E7.75': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.75_HC0.5_min400.h5ad"),
  'E8.0': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E8.0_HC0.5_min400.h5ad")
}

for stage,data in exp_data.items():
    data.obs.celltype = pd.Categorical(data.obs.celltype)

coord_data = {
  'E7.5': pd.read_csv('/rad/wuc/dash_data/spatial/spatial_coordinate.embryo_E7.5.csv', sep=' '),
  'E7.75': pd.read_csv('/rad/wuc/dash_data/spatial/spatial_coordinate.embryo_E7.75.csv', sep=' '),
  'E8.0': pd.read_csv('/rad/wuc/dash_data/spatial/spatial_coordinate.embryo_E8.0.csv', sep=' ')
}

for stage in exp_data.keys():
  exp_data[stage].obs[['x','y','z','x_flatten', 'y_flatten']] = coord_data[stage].loc[exp_data[stage].obs_names,['x','y','z','x_flatten', 'y_flatten']]

for stage in exp_data.keys():
  exp_data[stage].raw = exp_data[stage].copy()
  sc.pp.normalize_total(exp_data[stage], target_sum=1e4)
  sc.pp.log1p(exp_data[stage])
  exp_data[stage].obsm = {'X_spatial': coord_data[stage].loc[exp_data[stage].obs_names,['x','y','z']]}

for k,v in exp_data.items():
  v.obs.germ_layer = [i.replace('exe-ecto', 'ectoderm') for i in v.obs.germ_layer.values]

auc_data = {}
regulon_geneset = {}
for stage in ['E7.5', 'E7.75', 'E8.0']:
  h5 = h5py.File( '/rad/wuc/dash_data/spatial/matrix_data/%s_auc_mtx.h5' % (stage))
  auc_mtx = pd.DataFrame(h5['auc_matrix'], index=h5['cell'][:].astype(str), columns=h5['regulon'][:].astype(str))
  auc_data[stage] = sc.AnnData(X=auc_mtx, obs=exp_data[stage].obs)
  regulon_geneset[stage] = json.loads(h5['geneset'][()])
  h5.close()
del(auc_mtx)

genes_min_pval = pd.read_csv("/rad/wuc/dash_data/spatial/sparkX_res/genes_padj_minVal.csv",
                             header=[0], index_col=[0,1]) 


genes_all_pval = pd.read_csv("/rad/wuc/dash_data/spatial/sparkX_res/genes_padj_combine.csv",
                             header=[0,1], index_col=[0,1]) 
genes_all_pval = genes_all_pval.loc[:,[('ecto', 'adjustedPval'), ('meso', 'adjustedPval'), ('endo', 'adjustedPval'), ('all', 'adjustedPval')]]
genes_all_pval.columns = ['ecto p.adj', 'meso p.adj', 'endo p.adj', 'all p.adj']
genes_all_pval = genes_all_pval.groupby(level=0, group_keys=False
                                               ).apply(lambda x: x.sort_values(by='all p.adj'))
ctp_cmap = pd.read_csv("/rad/wuc/dash_data/spatial/celltype_cmap.csv")
ctp_cmap = dict(zip(ctp_cmap['celltype'], ctp_cmap['Epiblast']))

def show_feature_spatial_3D(adata, feature, embedding, cmap = None, sort=False, ascending=True, **kws):
    embedding = embedding.loc[adata.obs_names,:]
    if cmap is None:
        cmap = [(0.00, "rgb(244,244,244)"),
                (0.05, "rgb(244, 244, 244)"),
                (1.00, "rgb(34, 94, 168)")
                ]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['x', 'y', 'z'])
    pdf = pd.concat([pdf, adata[:,feature].to_df()], axis=1)
    if sort is True:
      pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter_3d(
        data_frame = pdf,
        x = 'x', y = 'y', z='z',
        color = feature,
        color_continuous_scale = cmap,
        **kws
      )
    plot.update_traces(marker_size=3,
                      marker_opacity = 0.8)
    plot.update_layout(
      margin=dict(l=10, r=10, t=30, b=0),
      uirevision='constant',
      coloraxis = {
        'colorbar' : {'tickformat': '4.2f'}
      }
    )
    
    plot.update_layout(
      scene = dict(
        xaxis = {
                'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        yaxis = {
                'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        zaxis = {
                'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        bgcolor = 'white'
      )
    )
    plot.add_trace(go.Mesh3d(x=[0,0], y=[0,0], z=[0,0]))
    return plot




# In[] js dict

dbc_css = "/home/wuc/dashapps/css/dbc.min.css"
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc_css],
          use_pages=False)


app.layout = dbc.Container(
  dbc.Col(
    [
      dcc.Store(id='store'),
      dcc.RangeSlider(min=-100, max=100, value=(-50,50), step=10, id='slider',
                      marks=None,tooltip={'placement': 'bottom', 'always_visible': True}),
      html.Div(id='div'),
    ]
  )
)

clientside_callback(
  '''
  function(value){
    var dict = {'min': value[0], 'max': value[1]}
    return dict
  }
  ''',
  Output('store', 'data'),
  Input('slider', 'value')
)

@app.callback(
  Output('div', 'children'),
  Input('store', 'data')
)
def show_store(data):
  return json.dumps(data)


if __name__ == "__main__":
  app.run(
  host='10.193.0.208',
  port='8051',
  # threaded=True,
  proxy=None,
  debug=True,
  use_reloader=False,
  # jupyter_mode='external',
)
# 

# In[] clientside_callback triggered_id
import dash_html_components as html
from dash import Dash
from dash.dependencies import Output, Input

# Create app.
app = Dash(prevent_initial_callbacks=True)
app.layout = html.Div([
    html.Button("Button 1", id="btn1"), html.Button("Button 2", id="btn2"), html.Div(id="log")
])
app.clientside_callback(
    """
    function(x, y){
        const triggered = dash_clientside.callback_context.triggered.map(t => t.prop_id);
        return "Hello from [" + triggered + "]";
    }
    """, Output("log", "children"), [Input("btn1", "n_clicks"), Input("btn2", "n_clicks")],
    prevent_initial_call=True)

if __name__ == '__main__':
    app.run_server(
      host='10.193.0.208',
      port='8051',
      debug=True
    )
# In[] plotly.js in clientside_callback

import dash_html_components as html
from dash import Dash
from dash.dependencies import Output, Input

# Create app.
app = Dash(prevent_initial_callbacks=True)
app.layout = html.Div([
  dbc.Button(id='button'),
  dcc.Graph(id='graph')
])
app.clientside_callback(
    """
    function(click){
      if(click){
        var trace1 = {
          x: [1, 2, 3, 4],
          y: [10, 15, 13, 17],
          mode: 'markers',
          type: 'scatter'
        }
        var trace2 = {
          x: [2, 3, 4, 5],
          y: [16, 5, 11, 9],
          mode: 'lines',
          type: 'scatter'
        }
        var trace3 = {
          x: [1, 2, 3, 4],
          y: [12, 9, 15, 12],
          mode: 'lines+markers',
          type: 'scatter'
        }
        var data = [trace1, trace2, trace3];
        return Plotly.newPlot(data)
      }
    }
    """,
    Output("graph", "figure"),
    Input("button", "n_clicks"),
    prevent_initial_call=True
  )

if __name__ == '__main__':
    app.run_server(
      host='10.193.0.208',
      port='8051',
      debug=True
    )


# In[] extendData api

import dash
import dash_html_components as html
import dash_core_components as dcc
import numpy as np

from dash.dependencies import Input, Output, State

adata = exp_data['E7.5']

dbc_css = "/home/wuc/dashapps/css/dbc.min.css"
app = Dash(
  __name__, 
  external_stylesheets=[dbc.themes.BOOTSTRAP, dbc_css],
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ]
)

app.layout = html.Div([
  dbc.Button(['Plot'],id='button'),
  dbc.Button(['JSON'], id='json'),
  dcc.Store(data={'obs': adata.obs[['x','y','z','germ_layer','celltype']].to_dict('index')}, id='JSONtoPlot'),
  dcc.Graph(figure={},id='graph'),
  dcc.Graph(figure={},id='graph2'),
  dcc.Store(data=ctp_cmap, id='ctp_cmap'),
  html.Div(id='text'),
])


clientside_callback(
  ClientsideFunction(
    namespace='clientside',
    function_name='test_expCtp_3Dfigure',
  ),
  Output("graph", "figure"),
  Output('graph2', 'figure'),
  Input("button", "n_clicks"),
  State('JSONtoPlot', 'data'),
  State('ctp_cmap', 'data'),
  prevent_initial_call=True,
)

@app.callback(
  Output('JSONtoPlot', 'data'),
  Input("json", "n_clicks"),
)
def generate_JSON(click):
  patch = Patch()
  # patch['obs'] = adata.obs[['x','y','z','germ_layer','celltype']].to_dict('index')
  singleExp = adata[:,'T'].X.toarray().reshape(-1)
  
  expFilter = adata.obs_names[singleExp >= 0]
  obsFilter = adata.obs_names[[ True if i in ['endoderm', 'mesoderm'] else False for i in adata.obs['germ_layer'] ]]
  patch['singleExp'] = singleExp
  patch['cellsExpFilter'] = expFilter
  patch['cellsObsFilter'] = obsFilter

  return patch


if __name__ == '__main__':
  app.run_server(
    host='10.193.0.208',
    port='8051',
    debug=True
  )
    
# %%
import dash
import dash_html_components as html
import dash_core_components as dcc
import numpy as np

from dash.dependencies import Input, Output, State

# Example data (a circle).
resolution = 1000
t = np.linspace(0, np.pi * 2, resolution)
x, y = np.cos(t), np.sin(t)
# Example app.
figure = dict(data=[{'x': [], 'y': []}], layout=dict(xaxis=dict(range=[-1, 1]), yaxis=dict(range=[-1, 1])))
app = dash.Dash(__name__, update_title=None)  # remove "Updating..." from title
app.layout = html.Div([dbc.Button(id="interval"),
    dcc.Graph(id='graph', figure=dict(figure)), 
    dcc.Store(id='offset', data=0), dcc.Store(id='store', data=dict(x=x, y=y, resolution=resolution)),
    html.Div(id='test')
])
app.clientside_callback(
    """
    function (n_intervals, data, offset) {
        offset = offset % data.x.length;
        const end = Math.min((offset + 10), data.x.length);
        return [[{x: [data.x.slice(offset, end)], y: [data.y.slice(offset, end)]}, [0], 100], end]
    }
    """,
    [Output('graph', 'extendData'), Output('offset', 'data')],
    [Input('interval', 'n_clicks')], [State('store', 'data'), State('offset', 'data')]
)

@app.callback(
  Output('test', 'children'),
  Input('graph', 'figure')
)
def test(fig):
  return 

if __name__ == '__main__':
  app.run_server(
    host='10.193.0.208',
    port='8051',
    debug=True
  )
# %%
adata.obs[['x','y','z','germ_layer','celltype']].to_dict('index')
# %%

adata = exp_data['E7.5']
feature = 'Cdx1'

dbc_css = "/home/wuc/dashapps/css/dbc.min.css"
app = Dash(
  __name__, 
  external_stylesheets=[dbc.themes.BOOTSTRAP, dbc_css],
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ]
)

range = {
        'x_min': -400, 'x_max': 480,
        'y_min': -410, 'y_max': 470,
        'z_min': 0, 'z_max': 700,
      }

app.layout = html.Div([
  dbc.Button(['Plot'],id='plot'),
  dbc.Button(['JSON'], id='json'),
  dcc.Store(data={'obs': adata.obs[['x','y','z','germ_layer','celltype']].to_dict('index')}, id='JSONtoPlot'),
  dcc.Graph(figure={},id='graph'),
  # dcc.Graph(figure={},id='graph2'),
  dcc.Store(data=ctp_cmap, id='ctp_cmap'),
  html.Div(id='text'),
])

@app.callback(
  Output('graph', 'figure'),
  Input('plot', 'n_clicks'),
)
def plot(click):
  pdf = pd.concat([adata[:,feature].to_df(), adata.obs.celltype], axis=1)
  counts = pdf.celltype.value_counts()
  counts = counts[counts>0]
  sorted_ctp = counts.index.to_list()
  pdf['celltype'] = pd.Categorical(pdf['celltype'].to_list(),
                                   categories=sorted_ctp[::-1])
  fig = px.violin(
    pdf, x=feature, y='celltype', color = 'celltype', 
    color_discrete_map=ctp_cmap, orientation='h', height=800,
  ).update_traces(
    side='positive', width=1.5, jitter=0.2, marker_size=2.5
  ).update_layout(
    plot_bgcolor = 'rgba(200,200,200,0.1)',
  ).update_yaxes(
    gridcolor='rgba(200,200,200,0.6)', gridwidth=1,
  ).update_xaxes(
    dtick=1, gridcolor='#ffffff', gridwidth=1, griddash='solid'
  )
  return fig
  

if __name__ == '__main__':
  app.run_server(
    host='10.193.0.208',
    port='8051',
    debug=True
  )
# %% collapse card

import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_antd_components as fac
import dash_bootstrap_components as dbc
import feffery_utils_components as fuc


app = Dash(
  __name__, 
  external_stylesheets=[dbc.themes.BOOTSTRAP],
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ]
)

images = [
  dmc.Image(radius="sm", src='https://images.unsplash.com/photo-1449824913935-59a10b8d2000?ixlib=rb-1.2.1&ixid=MnwxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8&auto=format&fit=crop&w=250&q=80'),
  dmc.Image(radius="sm", src='https://images.unsplash.com/photo-1444723121867-7a241cacace9?ixlib=rb-1.2.1&ixid=MnwxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8&auto=format&fit=crop&w=250&q=80'),
  dmc.Image(radius="sm", src='https://images.unsplash.com/photo-1444084316824-dc26d6657664?ixlib=rb-1.2.1&ixid=MnwxMjA3fDB8MHxwaG90by1wYWdlfHx8fGVufDB8fHx8&auto=format&fit=crop&w=250&q=80'),
]

card1 = fac.AntdCollapse(
    images,
    isOpen=False,
    title='Plot multi features',
    style={'width': 300, 'font-size': 16, 'font-family': 'Segoe UI',}
)

app.layout = html.Div([card1])


if __name__ == '__main__':
  app.run_server(
    host='10.193.0.208',
    port='8053',
    debug=True
  )

# %%
import time
import plotly.express as px
from dash_extensions.enrich import DashProxy, Output, Input, State, Serverside, html, dcc, \
    ServersideOutputTransform

app = DashProxy(transforms=[ServersideOutputTransform()])
app.layout = html.Div(
    [
        html.Button("Query data", id="btn"),
        dcc.Dropdown(id="dd"),
        dcc.Graph(id="graph"),
        dcc.Loading(dcc.Store(id="store"), fullscreen=True, type="dot"),
    ]
)

@app.callback(Output("store", "data"), Input("btn", "n_clicks"), prevent_initial_call=True)
def query_data(n_clicks):
    time.sleep(3)  # emulate slow database operation
    return Serverside(px.data.gapminder())  # no JSON serialization here

@app.callback(Output("dd", "options"),  Output("dd", "value"), Input("store", "data"), prevent_initial_call=True)
def update_dd(df):
    options = [{"label": column, "value": column} for column in df["year"]]   # no JSON de-serialization here
    return options, options[0]['value']

@app.callback(Output("graph", "figure"), [Input("dd", "value"), State("store", "data")], prevent_initial_call=True)
def update_graph(value, df):
    df = df.query("year == {}".format(value))  # no JSON de-serialization here
    return px.sunburst(df, path=["continent", "country"], values="pop", color="lifeExp", hover_data=["iso_alpha"])

if __name__ == "__main__":
    app.run_server(
      host = '10.193.0.208',
      port = '8053',
      debug = True
    )
# %%
