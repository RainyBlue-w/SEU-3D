## In[] env

from math import isnan
import math
from dash import Dash, dcc, html, dash_table, Input, Output, callback, no_update, State, Patch, DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash.exceptions import PreventUpdate
import dash_daq as daq
import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_antd_components as fac
import dash_bootstrap_components as dbc
import feffery_utils_components as fuc
from dash_extensions.enrich import Output, Input, html, DashProxy, LogTransform, DashLogger

import plotly.express as px
import plotly.graph_objects as go
import plotly
from plotly.subplots import make_subplots
from plotnine import *
import plotnine.options

from PIL import Image
import scanpy as sc
import os
import pandas as pd
import dask.dataframe as dd
import numpy as np
# import loompy as lp
import h5py
import json
import time

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


## In[] data

exp_data = {
  'E7.5': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.5_HC0.5_min400.h5ad"),
  'E7.75': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.75_HC0.5_min400.h5ad"),
  'E8.0': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E8.0_HC0.5_min400.h5ad")
}

# for stage, data in exp_data.items():
  

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

for stage in exp_data.keys():
  exp_data[stage]= exp_data[stage][exp_data[stage].obs_names.sort_values()]

for k,v in exp_data.items():
  v.obs.germ_layer = [i.replace('exe-ecto', 'ectoderm') for i in v.obs.germ_layer.values]

auc_data = {}
regulon_geneset = {}
for stage in ['E7.5', 'E7.75', 'E8.0']:
  h5 = h5py.File( '/rad/wuc/dash_data/spatial/matrix_data/%s_auc_mtx.h5' % (stage))
  auc_mtx = pd.DataFrame(h5['auc_matrix'], index=h5['cell'][:].astype(str), columns=h5['regulon'][:].astype(str))
  auc_data[stage] = sc.AnnData(X=auc_mtx.loc[exp_data[stage].obs_names,:], obs=exp_data[stage].obs)
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


# In[] functions:
def show_expViolin(adata, feature, **kws):
  data = adata[:,feature].to_df()[feature]
  # data = data[data>0]
  fig = go.Figure(
    data = go.Violin(
      x=data, y0=f'{feature}({len(data)})', box_visible=True, 
      line_color='black', meanline_visible=True,
      fillcolor='lightseagreen', opacity=0.6,
      pointpos = 1.5, jitter=0.2, **kws
    )
  )
  
  fig.update_traces(orientation='h', side='positive', points=False, marker_size=2.5)
  fig.update_layout(xaxis_showgrid=False, yaxis_showgrid=True)
  return fig

def show_ctpExpViolin(adata, feature, **kws):
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

def show_multiFeatures_expViolin(adata, features_dict, color_dict=None,**kws):
  if not color_dict:
    color_dict = {'R': 'tomato', 'G': 'springgreen', 'B': 'skyblue'}
  fig = go.Figure()
  
  tmp = {}
  for key,value  in features_dict.items():
      if value:
          tmp[key] = value
  features_dict = tmp

  for color,feature in features_dict.items():
    data = adata[:,feature].to_df()[feature]
    # data = data[data>0]
    color = color_dict[color]
    fig.add_trace(
      go.Violin(
        x=data, y0=f'{feature}({len(data)})', box_visible=True, 
        line_color='black', meanline_visible=True,
        fillcolor=color, opacity=0.6,
        pointpos = 1.5, jitter=0.2, **kws
      )
    )
  fig.update_traces(orientation='h', side='positive', points=False, 
                    width=1.2)
  fig.update_layout(showlegend=False,xaxis_showgrid=False, yaxis_showgrid=True)
  return fig

def show_multiFeatures_ctpExpViolin(adata, features_dict, *kws):
  
  from plotly.subplots import make_subplots
  
  tmp = {}
  for key,value  in features_dict.items():
      if value:
          tmp[key] = value
  features_dict = tmp
  features = list(features_dict.values())
  ngenes = len(features_dict.keys())
  fig = make_subplots(1, ngenes)

  pdf = pd.concat([adata[:,features].to_df(), adata.obs.celltype], axis=1)
  pdf = pdf.melt(id_vars='celltype')
  pdf = pdf.rename(columns = {'variable': 'Gene', 'value': 'expression'})
  # pdf = pdf[pdf['expression']>0]

  pdf.celltype = pd.Categorical(pdf.celltype, ordered=True)
  # counts = pdf.groupby('Gene').apply(lambda x: x.value_counts())

  fig = px.violin(
    pdf, x='expression', y='celltype', color = 'celltype', 
    color_discrete_map=ctp_cmap, orientation='h', height=800,
    animation_frame='Gene', points='outliers'
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

def vector_to_rgba(v):
  color = list(v.keys())
  color = [str(math.ceil(v[i])) if i in color else '244' for i in ['R', 'G', 'B'] ]
  if(all([ i=='244' for i in color])):
    rgba = 'rgba(244,244,244,1)'
  else:
    rgba = f'rgba({color[0]}, {color[1]}, {color[2]}, 1)'
    
  return rgba

def multiGenes_show_color(adata, genes_dict):
  import numpy
  tmp = {}
  for key,value  in genes_dict.items():
      if value:
          tmp[key] = value
  genes_dict = tmp
  colors = list(genes_dict.keys())
  others = [i for i in ['R', 'G', 'B'] if i not in colors]
  genes = list(genes_dict.values())

  exp = adata[:, genes].to_df()
  exp.columns = colors

  delta = 244 - exp.div(exp.max(axis=0), axis=1)*244
  delta[others] = 244

  def color_geoMean(a,b):
    a = 244-a
    b = 244-b
    geoMean = numpy.sqrt((a**2+b**2)/2)
    # geoMean = ((a**3+b**3)/2)**(1/3)
    color = 244 - geoMean
    return color
  def mean(a,b, c=None):
    if c:
      return (a+b+c)/3
    else:
      return (a+b)/2

  if len(colors)==1:
    color = pd.DataFrame({
        colors[0] : 244,
        others[0] : delta[colors[0]],
        others[1] : delta[colors[0]],
    })
  elif len(colors)==2:
    color = pd.DataFrame({
        colors[0] : delta[colors[1]],
        colors[1] : delta[colors[0]],
        others[0] : color_geoMean(delta[colors[1]],delta[colors[0]]),
    })
  elif len(colors)==3:
    color = pd.DataFrame({
        'R' : color_geoMean(delta['G'], delta['B']),
        'G' : color_geoMean(delta['R'], delta['B']),
        'B' : color_geoMean(delta['R'], delta['G']),
    })
  
  color['RGBA'] = color.apply(vector_to_rgba, axis=1)
  return color['RGBA']

def cal_moran_3D(adata):
  tmp = adata.copy()
  sq.gr.spatial_neighbors(tmp, spatial_key='X_spatial')
  sq.gr.spatial_autocorr(tmp, mode='moran', n_jobs=1)
  df = tmp.uns['moranI'][['I']]
  df.columns = ["Moran's I"]
  return df

# In[] app

app = DashProxy(
  __name__, 
  external_stylesheets=[
    dbc.themes.BOOTSTRAP
  ],
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ],
  transforms=[
    LogTransform()
  ]
)


header = dbc.NavbarSimple(
    [
        dbc.DropdownMenu(
            children=[
                # dbc.DropdownMenuItem('spatial', href='/'),
                # dbc.DropdownMenuItem("atlas", href='/atlas'),
                # dbc.DropdownMenuItem("reik", href='/reik'),
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

# In[] widgets

SET_STORE_JSONtoPlot_3D = html.Div(
  [
    dcc.Store(data={}, id='STORE_obs_3D'),
    dcc.Store(data={}, id='STORE_cellsObsFilter_3D'),
    dcc.Store(data={}, id='STORE_cellsExpFilter_3D'),
    dcc.Store(data={}, id='STORE_singleExp_3D'),
    dcc.Store(data={}, id='STORE_multiExp_3D'),
    dcc.Store(data={}, id='STORE_mixedColor_3D'),
    dcc.Store(data=False, id='STORE_ifmulti_3D'),
    dcc.Store(data=ctp_cmap, id='STORE_ctpCmap_3D'),
    dcc.Store(id='STORE_ctpToShow_3D'),
    dcc.Store(id='STORE_cellsCtpFilter_3D'),
    dcc.Store(id='test'),
  ]
)

SET_STORE_Ranges_3D = html.Div(
  [
    dcc.Store(
      data= {
        'x_min': -400, 'x_max': 480,
        'y_min': -410, 'y_max': 470,
        'z_min': 0, 'z_max': 700,
      },
      id='STORE_previewRange_3D'
    ),
    dcc.Store(id='STORE_sliceRange_3D'),
    dcc.Store(
      data= {
        'x_min': -400, 'x_max': 480,
        'y_min': -410, 'y_max': 470,
        'z_min': 0, 'z_max': 700,
      },
      id='STORE_maxRange_3D'),
  ]
)

CHECKLIST_germLayer_3D = html.Div(
    [
        # dbc.Label("Germ layers to show:"),
        dbc.Checklist(
            options=[
                {"label": " Ectoderm", "value": 'ectoderm'},
                {"label": " Mesoderm", "value": 'mesoderm'},
                {"label": " Endoderm", "value": 'endoderm'},
            ],
            value=['ectoderm', 'mesoderm', 'endoderm'],
            id="CHECKLIST_germLayer_3D",
            switch=True
        ),
    ],
    className="mb-4",
)

SWITCH_hideZero_3D = html.Div([
  dbc.Row(
    [
        dbc.Col(
            dbc.Label("Hide zero:"),
            width=6
            ),
        dbc.Col(
            daq.BooleanSwitch(id='SWITCH_hideZero_3D', on=False),
            width=4
            )
    ]
  )
])

SWITCH_hideAxes_3D = html.Div([
  dbc.Row(
    [
      dbc.Col(
        dbc.Label('Hide axes:'),
        width=6,
      ),
      dbc.Col(
        daq.BooleanSwitch(on=False, id='SWITCH_hideAxes_3D'),
        width=4
      )
    ]
  )
])

SET_topViewer_controler_3D = html.Div([
  dbc.Row(
    [
      dbc.Col([
        dbc.Label('x range'),
        dcc.RangeSlider(-400,480,10, value=(-400,480), id='SLIDER_Xrange_3D',
                        marks=None,tooltip={'placement': 'bottom', 'always_visible': True})
      ], width=2),
      dbc.Col([
        dbc.Label('y range'),
        dcc.RangeSlider(-410, 470,10, value=(-410,470), id='SLIDER_Yrange_3D',
                        marks=None,tooltip={'placement': 'bottom', 'always_visible': True})
      ], width=2),
      dbc.Col([
        dbc.Label('z range'),
        dcc.RangeSlider(0, 700,10, value=(0,700), id='SLIDER_Zrange_3D',
                        marks=None,tooltip={'placement': 'bottom', 'always_visible': True})
      ], width=2),
      dbc.Col([
        dbc.Row([
          dbc.Col(dbc.Label('Preview'), width=4),
          dbc.Col(daq.BooleanSwitch(on=False, id='SWITCH_previewBox_3D'),width=4)
        ]),
        dbc.Row([
          dbc.Col(dbc.Button('Slice', color='danger', id='BUTTON_slice_3D'),width=4),
          dbc.Col(dbc.Button('Recover', color='success', id='BUTTON_recover_3D'),width=4),
        ])
      ], width=2),
      dbc.Col([
        CHECKLIST_germLayer_3D
      ], width=2),
      dbc.Col([
        SWITCH_hideZero_3D,
        SWITCH_hideAxes_3D
      ], width=2)
    ],
    style = {'height': '8vh'}
  ),
  SET_STORE_Ranges_3D,
], className='mb-4')


# In[] tab

spatial_tab_plotFeature3D = dbc.Tab(
  [dbc.Row([
    # options
    dbc.Col([
      html.Div([
        # Basic options
        fac.AntdCollapse(
          [
            html.Div([
              dbc.Row([
                dbc.Col([
                  dbc.Label("Feature type"),
                  dcc.Dropdown(
                      ['Gene', 'Regulon'],
                      'Gene',
                      id="DROPDOWN_featureType_3D",
                      clearable=False,
                      searchable=True,
                  ),
                ], width=6),
                dbc.Col([
                  dbc.Label("Stage"),
                  dcc.Dropdown(
                      ['E7.5', 'E7.75', 'E8.0'],
                      'E7.75',
                      id="DROPDOWN_stage_3D",
                      clearable=False,
                      searchable=True,
                  ),
                ], width=6),
              ]),
              
            ]),
            
          ],
          isOpen=True,
          title='Select data',
          forceRender=True,
          style={'font-size': 16, 'font-family': 'Segoe UI',}
        ),
        # Slicer
        # fac.AntdCollapse(
          
        # ),
        # Single gene
        fac.AntdCollapse(
          [
            html.Div([
              dbc.Row(
                [
                  dbc.Col([
                    dcc.Dropdown(
                      options = exp_data['E7.75'].var_names,
                      value = 'Cdx1',
                      id="DROPDOWN_singleName_3D",
                      clearable=False
                    ),
                  ], width=8),
                  dbc.Col([
                    dbc.Button('Plot', id='BUTTON_singlePlot_3D', n_clicks=0, color='dark', disabled=False),
                  ], width=4),
                ],className="mb-4",
              )
            ])
          ],
          isOpen=True,
          title='Plot single feature',
          forceRender=True,
          style={'font-size': 16, 'font-family': 'Segoe UI',}
        ),
        # Multi genes
        fac.AntdCollapse(
          [
            html.Div([
              dbc.Col([
                dbc.Row([
                  dbc.Col(dcc.Dropdown(options = exp_data['E7.75'].var_names,
                                      id='DROPDOWN_multiNameR_3D'),
                          width=10),
                  dbc.Col(dbc.Badge("R", color='danger'),width=2),
                ]),
                dbc.Row([
                  dbc.Col(dcc.Dropdown(options = exp_data['E7.75'].var_names,
                                      id='DROPDOWN_multiNameG_3D'), 
                          width=10),
                  dbc.Col(dbc.Badge("G", color='success'),width=2),
                ]),
                dbc.Row([
                  dbc.Col(dcc.Dropdown(options = exp_data['E7.75'].var_names,
                                      id='DROPDOWN_multiNameB_3D'), 
                          width=10),
                  dbc.Col(dbc.Badge("B", color='primary'),width=2),
                ]),
                dbc.Row([
                  dbc.Button('Plot', id='BUTTON_multiPlot_3D', n_clicks=0, color='dark', disabled=False),
                ], justify='center'),
                dcc.Store(id='STORE_multiNameInfo_3D')
              ],)
            ])
          ],
          isOpen=False,
          title='Plot multi features',
          forceRender=True,
          style={'font-size': 16, 'font-family': 'Segoe UI',}
        ),
        # Moran
        fac.AntdCollapse(
          [
            html.Div([
              dbc.Row([
                dbc.Col(dbc.Button('calculate\nSVGs', id='BUTTON_calMoran_3D', color='dark', disabled=False), width=6),
                dbc.Col(dbc.Button('Show\nresult', id='BUTTON_showMoran_3D', color='primary', disabled=False), width=6),
              ]),
              dbc.Offcanvas(
                [dash_table.DataTable(
                  id='DATATABLE_moranRes_3D',
                  sort_action="native", page_action='native', filter_action="native",
                  page_current= 0, page_size= 20, fill_width=True,
                  style_cell={'textAlign': 'center'},
                  style_table={'overflowX': 'auto'},
                  # style_cell={'padding-right': '10px', 'padding-left': '10px',
                  # 'text-align': 'center', 'marginLeft': 'auto', 'marginRight': 'auto'})],
                )],
                title = 'SVGs:',
                placement='end', scrollable=True, backdrop=False, is_open=False,
                id = 'OFFCANVAS_moranRes_3D',
              ),
            ]),
          ],
          isOpen=False,
          title='Calculate SVG(moran)',
          style={'font-size': 16, 'font-family': 'Segoe UI',}
        )
      ], style = {'position':'fixed', 'width':'30vh', 'overflowY': 'overlay', 'maxHeight': '90vh'}),
    ], width=2),
    # viewer
    dbc.Col([
      SET_topViewer_controler_3D,
      SET_STORE_JSONtoPlot_3D,
      # scatter3d
      dbc.Row([
        dbc.Col([
          dcc.Graph(figure={}, id="FIGURE_3Dexpression", 
                    style={'height': "80vh"}, ),
        ],align = "center", width=5),
        dbc.Col([
          dcc.Graph(figure={}, id="FIGURE_3Dcelltype", 
                    style={'height': "80vh"}),
        ],align = "center", width=7)
      ]),
      # violin
      dbc.Row([
        dbc.Col([
          dbc.Label( 'Normalized expression in all celltypes(left)'),
          dbc.Label('and in each celltype(right):'),
          dmc.LoadingOverlay(dcc.Graph(figure={}, id="FIGURE_expViolin_3D"))
        ], align='center', width=4),
        dbc.Col([
          dmc.LoadingOverlay(dcc.Graph(figure={}, id="FIGURE_ctpViolin_3D"))
        ], align='center', width=8)
      ])
    ],width=10),
  ])],
  label = "Plot feature(3D)",
  tab_id = "spatial_tab_plotFeature3D",
)

spatial_tabs = dbc.Card(
    dbc.Tabs(
        [spatial_tab_plotFeature3D],
        active_tab = "spatial_tab_plotFeature3D",  
        id = "spatial_tabs",
    ),
)

# In[] callbacks

# update_nameOptions
@app.callback(
  Output('DROPDOWN_singleName_3D', 'options'),
  Input('DROPDOWN_singleName_3D', 'search_value'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value')
)
def update_nameOptions_single_3D(search, featureType, stage):
  if not search:
    raise PreventUpdate
  
  if featureType == 'Gene':
    if not search:
      opts = exp_data[stage].var_names
    else:
      opts = exp_data[stage].var_names[exp_data[stage].var_names.str.startswith(search)].sort_values()
  elif featureType == 'Regulon':
    if not search:
      opts = auc_data[stage].var_names
    else:
      opts = auc_data[stage].var_names[auc_data[stage].var_names.str.startswith(search)].sort_values()
  
  return opts

@app.callback(
  Output('DROPDOWN_multiNameR_3D', 'options'),
  Input('DROPDOWN_multiNameR_3D', 'search_value'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value')
)
def update_nameOptions_multiR_3D(search, featureType, stage):
  if not search:
    raise PreventUpdate
  
  if featureType == 'Gene':
    if not search:
      opts = exp_data[stage].var_names
    else:
      opts = exp_data[stage].var_names[exp_data[stage].var_names.str.startswith(search)].sort_values()
  elif featureType == 'Regulon':
    if not search:
      opts = auc_data[stage].var_names
    else:
      opts = auc_data[stage].var_names[auc_data[stage].var_names.str.startswith(search)].sort_values()
  
  return opts

@app.callback(
  Output('DROPDOWN_multiNameG_3D', 'options'),
  Input('DROPDOWN_multiNameG_3D', 'search_value'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value')
)
def update_nameOptions_multiG_3D(search, featureType, stage):
  if not search:
    raise PreventUpdate
  
  if featureType == 'Gene':
    if not search:
      opts = exp_data[stage].var_names
    else:
      opts = exp_data[stage].var_names[exp_data[stage].var_names.str.startswith(search)].sort_values()
  elif featureType == 'Regulon':
    if not search:
      opts = auc_data[stage].var_names
    else:
      opts = auc_data[stage].var_names[auc_data[stage].var_names.str.startswith(search)].sort_values()
  
  return opts

@app.callback(
  Output('DROPDOWN_multiNameB_3D', 'options'),
  Input('DROPDOWN_multiNameB_3D', 'search_value'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value')
)
def update_nameOptions_multiB_3D(search, featureType, stage):
  if not search:
    raise PreventUpdate
  
  if featureType == 'Gene':
    if not search:
      opts = exp_data[stage].var_names
    else:
      opts = exp_data[stage].var_names[exp_data[stage].var_names.str.startswith(search)].sort_values()
  elif featureType == 'Regulon':
    if not search:
      opts = auc_data[stage].var_names
    else:
      opts = auc_data[stage].var_names[auc_data[stage].var_names.str.startswith(search)].sort_values()
  
  return opts

# store_multiNameInfo
app.clientside_callback(
  '''
  function(R,G,B){
    var dict = {'R': R, 'G': G, 'B': B};
    return dict
  }
  ''',
  Output('STORE_multiNameInfo_3D', 'data'),
  Input('DROPDOWN_multiNameR_3D', 'value'),
  Input('DROPDOWN_multiNameG_3D', 'value'),
  Input('DROPDOWN_multiNameB_3D', 'value')
)

# store_previewRange
app.clientside_callback(
  ClientsideFunction(
    namespace='plotFunc_3Dtab',
    function_name='store_previewRange',
  ),
  Output('STORE_previewRange_3D', 'data'),
  Input('SLIDER_Xrange_3D', 'value'),
  Input('SLIDER_Yrange_3D', 'value'),
  Input('SLIDER_Zrange_3D', 'value'),
)

# store_sliceRange
app.clientside_callback(
  '''
  function(slice, recover, previewRange, maxRange){
    const id = dash_clientside.callback_context.triggered.map(t => t.prop_id)
    console.log(id)
    if(id.includes('BUTTON_slice_3D.n_clicks')){
      return previewRange
    } else {
      return maxRange
    }
  }
  ''',
  Output('STORE_sliceRange_3D', 'data'),
  Input('BUTTON_slice_3D', 'n_clicks'),
  Input('BUTTON_recover_3D', 'n_clicks'),
  State('STORE_previewRange_3D', 'data'),
  State('STORE_maxRange_3D', 'data'),
)

# store_cellsInfo_forJSONtoPlot (download: ~2.5M,500ms ; compute 250ms)
@app.callback(
  Output('STORE_cellsObsFilter_3D', 'data'),
  
  Input('STORE_sliceRange_3D', 'data'),
  Input('CHECKLIST_germLayer_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value'),
  Input('DROPDOWN_featureType_3D', 'value'),
)
def store_cellsInfo_forJSONtoPlot_3D(sliceRange, germs, stage, featureType):

  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_data[stage]
  else:
    raise PreventUpdate
  
  obs = adata.obs[['x','y','z','germ_layer','celltype']]

  if_inSliceRange = ( 
                      (obs['x'] <= sliceRange['x_max']) & 
                      (obs['x'] >= sliceRange['x_min']) & 
                      (obs['y'] <= sliceRange['y_max']) & 
                      (obs['y'] >= sliceRange['y_min']) & 
                      (obs['z'] <= sliceRange['z_max']) & 
                      (obs['z'] >= sliceRange['z_min'])
                    )
  
  if len(germs) >= 1:
    if_ingerms = [ True if i in germs else False for i in obs['germ_layer'] ]
  else:
    raise PreventUpdate

  obsnames_filt = adata.obs_names[if_inSliceRange & if_ingerms]
  return obsnames_filt
  
@app.callback(
  Output('STORE_obs_3D', 'data'),
  Input('DROPDOWN_stage_3D', 'value'),
  Input('DROPDOWN_featureType_3D', 'value'),
)
def store_cellsObs_forJSONtoPlot_3D(stage, featureType):
  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_data[stage]
  else:
    raise PreventUpdate
  
  obs = adata.obs[['x','y','z','germ_layer','celltype']]
  return obs.to_dict('index')

# store_expInfo_forJSONtoPlot (download: <0.43M,<80ms; compute 320ms)
@app.callback(
  Output('STORE_cellsExpFilter_3D', 'data'),
  Output('STORE_singleExp_3D', 'data'),
  Output('STORE_multiExp_3D', 'data'),
  Output('STORE_ifmulti_3D', 'data'),
  Output('STORE_mixedColor_3D', 'data'),
  
  Input('BUTTON_singlePlot_3D', 'n_clicks'),
  Input('BUTTON_multiPlot_3D', 'n_clicks'),
  Input('SWITCH_hideZero_3D', 'on'),
  Input('DROPDOWN_stage_3D', 'value'),
  Input('DROPDOWN_featureType_3D', 'value'),
  
  State('DROPDOWN_singleName_3D', 'value'),
  State('STORE_multiNameInfo_3D', 'data'),
  State('STORE_ifmulti_3D', 'data'),

)
def store_expInfo_forJSONtoPlot_3D(sclick, mclick, hideZero, stage, featureType, sname, minfo, ifmulti):

  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_data[stage]
  else:
    raise PreventUpdate

  
  def return_single():
    ifmulti = False
    exp = adata[:,sname].to_df()
    if hideZero:
      cellsExpFilter = adata.obs_names[(exp>0)[sname]].to_list()
    else:
      cellsExpFilter = adata.obs_names.to_list()
    return (ifmulti, exp, cellsExpFilter)
  
  def return_multi():
    ifmulti = True
    mixColor = dict(zip( adata.obs_names, multiGenes_show_color(adata, minfo)))
    if hideZero:
      cellsExpFilter = adata.obs_names[[ i != 'rgba(244,244,244,1)' for i in mixColor.values()]].to_list()
    else:
      cellsExpFilter = adata.obs_names.to_list()
    return (ifmulti, [], cellsExpFilter, mixColor) 
  
  def return_multiExp():
    tmp = {}
    for key,value in minfo.items():
      if value:
        tmp[key] = value
    colors = list(tmp.keys())
    genes = list(tmp.values())

    exp = adata[:, genes].to_df()
    exp.columns = colors
    exp.to_dict('index')

    return exp
  
  btn_id = ctx.triggered_id
  if btn_id:
    if 'DROPDOWN_stage_3D' in btn_id or 'DROPDOWN_featureType_3D' in btn_id:
      if not ifmulti:
        ifmulti,exp,cellsExpFilter = return_single()
        exp = exp.to_dict('index')
        return (cellsExpFilter, exp, no_update, ifmulti, no_update)
      else:
        ifmulti,_,cellsExpFilter,mixcolor = return_multi()
        return (cellsExpFilter, no_update, no_update, ifmulti, mixcolor)

    elif 'BUTTON_singlePlot_3D' in btn_id:
      ifmulti,exp,cellsExpFilter = return_single()
      exp = exp.to_dict('index')
      if hideZero:
        return (cellsExpFilter, exp, no_update, ifmulti, no_update)
      else:
        return (no_update, exp, no_update, ifmulti, no_update)
    
    elif 'BUTTON_multiPlot_3D' in btn_id:
      ifmulti,_,cellsExpFilter,mixcolor = return_multi()
      # exp = return_multiExp()
      if hideZero:
        return (cellsExpFilter, no_update, no_update, ifmulti, mixcolor)
      else:
        return (no_update, no_update, no_update, ifmulti, mixcolor)
    
    elif 'SWITCH_hideZero_3D' in btn_id:
      
      if not hideZero:
        cellsExpFilter = adata.obs_names.to_list()
        return (cellsExpFilter, no_update, no_update, no_update, no_update)
      
      else:
        if not ifmulti:
          _,_,cellsExpFilter = return_single()
          return (cellsExpFilter, no_update, no_update, no_update, no_update)
        else:
          _,_,cellsExpFilter,_ = return_multi()
          return (cellsExpFilter, no_update, no_update, no_update, no_update)

  else:
      ifmulti,exp,cellsExpFilter = return_single()
      exp = exp.to_dict('index')
      return (cellsExpFilter, exp, no_update, ifmulti, no_update)

# plot_3Dfigure_exp
app.clientside_callback(
  ClientsideFunction(
    namespace='plotFunc_3Dtab',
    function_name='exp_3Dscatter',
  ),
  Output("FIGURE_3Dexpression", "figure"),
  Input('STORE_obs_3D', 'data'),
  Input('STORE_cellsObsFilter_3D', 'data'),
  Input('STORE_cellsExpFilter_3D', 'data'),
  Input('STORE_cellsCtpFilter_3D', 'data'),
  Input('STORE_singleExp_3D', 'data'),
  Input('STORE_ifmulti_3D', 'data'),
  Input('STORE_mixedColor_3D', 'data'),
  State('SWITCH_hideAxes_3D', 'on'),
  State('SWITCH_previewBox_3D', 'on'),
  State('STORE_previewRange_3D', 'data'),
)

# plot_3Dfigure_ctp

app.clientside_callback(
  ClientsideFunction(
    namespace='plotFunc_3Dtab',
    function_name='ctp_3Dscatter',
  ),
  Output("FIGURE_3Dcelltype", "figure"),
  Input('STORE_obs_3D', 'data'),
  Input('STORE_cellsObsFilter_3D', 'data'),
  Input('STORE_cellsExpFilter_3D', 'data'),
  State('SWITCH_hideAxes_3D', 'on'),
  State('SWITCH_previewBox_3D', 'on'),
  State('STORE_previewRange_3D', 'data'),
  State('STORE_ctpCmap_3D', 'data'),
)

# sync layout between exp and ctp figure
@app.callback(
  Output("FIGURE_3Dexpression", "figure", allow_duplicate=True),
  Output("FIGURE_3Dcelltype", "figure", allow_duplicate=True),
  Input("FIGURE_3Dexpression", "relayoutData"),
  Input("FIGURE_3Dcelltype", "relayoutData"),
  prevent_initial_call=True
)
def update_relayout(expLayout, ctpLayout):
  tid = ctx.triggered_id
  patch = Patch()
  if tid == 'FIGURE_3Dexpression':
    patch['layout']['scene']['camera'] = expLayout['scene.camera']
    return no_update, patch
  if tid == 'FIGURE_3Dcelltype':
    patch['layout']['scene']['camera'] = ctpLayout['scene.camera']
    return patch, no_update

# sync restyle between exp and ctp figure
@app.callback(
  Output('STORE_ctpToShow_3D', 'data'),
  Output('STORE_cellsCtpFilter_3D', 'data'),
  
  Input('DROPDOWN_stage_3D', 'value'),
  Input('BUTTON_singlePlot_3D', 'n_clicks'),
  Input('BUTTON_multiPlot_3D', 'n_clicks'),
  Input('STORE_cellsObsFilter_3D', 'data'),
  Input('STORE_cellsExpFilter_3D', 'data'),
  State('STORE_ctpToShow_3D', 'data'),
  Input('SWITCH_hideZero_3D', 'on'),
  Input('FIGURE_3Dcelltype', 'restyleData'),
)
def sync_restyle_3D(stage, splot, mplot, cellsObsFilter, cellsExpFilter, ctpNow, hideZero, restyle):
  
  cellsfilter = list(set(cellsExpFilter) & set(cellsObsFilter))
  cellsfilter.sort()
  adata =  exp_data[stage][cellsfilter]
  celltype_all = list(adata.obs['celltype'].unique())
  
  id = ctx.triggered_id
 
# --- for initial ---
  if (not ctpNow) or ('STORE_cellsObsFilter_3D' in id) or ('SWITCH_hideZero_3D' in id):
    ctpNow = celltype_all
    restyle = [{'visible': [True]}, [0]]

    
  if (('BUTTON_singlePlot_3D' in id) or ('BUTTON_multiPlot_3D' in id)) and hideZero:
    ctpNow = celltype_all
    restyle = [{'visible': [True]}, [0]]

    
# ------------------
  if(restyle):

    for index,order in enumerate(restyle[1]):
      if index != len(celltype_all):
        ctpNow[order] = None if restyle[0]['visible'][index] == 'legendonly' else celltype_all[order]


  cells = adata.obs_names[[True if i in ctpNow else False for i in adata.obs.celltype ]].to_list()
  
  return (ctpNow, cells)

# hide axes
@app.callback(
  Output("FIGURE_3Dexpression", "figure", allow_duplicate=True),
  Output("FIGURE_3Dcelltype", "figure", allow_duplicate=True),
  Input('SWITCH_hideAxes_3D', 'on'),
  prevent_initial_call=True
)
def hideAxes_3D(hideAxes):
  patch = Patch()
  if hideAxes:
    patch['layout']['scene']['xaxis']['visible'] = False
    patch['layout']['scene']['yaxis']['visible'] = False
    patch['layout']['scene']['zaxis']['visible'] = False
  else: 
    patch['layout']['scene']['xaxis']['visible'] = True
    patch['layout']['scene']['yaxis']['visible'] = True
    patch['layout']['scene']['zaxis']['visible'] = True
  return patch, patch

# show preview Box
@app.callback(
  Output('FIGURE_3Dexpression', 'figure', allow_duplicate=True),
  Output('FIGURE_3Dcelltype', 'figure', allow_duplicate=True),
  Input('SWITCH_previewBox_3D', 'on'),
  Input('STORE_previewRange_3D', 'data'),
  prevent_initial_call=True,
)
def update_previewBox(showBox, preRange):
  patch = Patch()
  if showBox:
    patch['data'][-1] = {
                    'x': [preRange['x_min'], preRange['x_min'], preRange['x_min'], preRange['x_min'],
                          preRange['x_max'], preRange['x_max'], preRange['x_max'], preRange['x_max']],
                    'y': [preRange['y_min'], preRange['y_max'], preRange['y_min'], preRange['y_max'],
                          preRange['y_min'], preRange['y_max'], preRange['y_min'], preRange['y_max']],
                    'z': [preRange['z_min'], preRange['z_min'], preRange['z_max'], preRange['z_max'],
                          preRange['z_min'], preRange['z_min'], preRange['z_max'], preRange['z_max']],
                    'i': [0, 1, 0, 0, 0, 0, 2, 2, 7, 7, 7, 7],
                    'j': [1, 2, 4, 1, 4, 2, 3, 6, 4, 4, 1, 1],
                    'k': [2, 3, 5, 5, 6, 6, 7, 7, 6, 5, 3, 5],
                    'color': 'black', 'opacity': 0.60, 'type': 'mesh3d'
                  }
  else:
    patch['data'][-1] = {
                    'x': [], 'y': [], 'z': [], 'i': [], 'j': [], 'k': [],
                    'type': 'mesh3d', 'color': 'black', 'opacity': 0.60
                  }

  return patch, patch

# violin plot
@app.callback(
  Output('FIGURE_expViolin_3D', 'figure'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value'),
  Input('STORE_cellsExpFilter_3D', 'data'),
  Input('STORE_cellsObsFilter_3D', 'data'),
  Input('STORE_cellsCtpFilter_3D', 'data'),
  Input('STORE_ifmulti_3D', 'data'),
  Input('BUTTON_singlePlot_3D', 'n_clicks'),
  Input('BUTTON_multiPlot_3D', 'n_clicks'),
  State('DROPDOWN_singleName_3D', 'value'),
  State('STORE_multiNameInfo_3D', 'data'),
  backgroud = True,
  manager = background_callback_manager
)
def update_spatial_plotFeature3D_expViolin(featureType, stage, cellsExp, cellsObs, cellsCtp, ifmulti, splot, mplot, sname, minfo):
  
  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]
  
  cells_to_use = list(set(cellsObs) & set(cellsExp) & set(cellsCtp))
  adata = adata[cells_to_use]

  if not ifmulti:
    fig = show_expViolin(adata, sname)
  else:
    fig = show_multiFeatures_expViolin(adata, minfo)

  return fig

@app.callback(
  Output('FIGURE_ctpViolin_3D', 'figure'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value'),
  Input('STORE_cellsExpFilter_3D', 'data'),
  Input('STORE_cellsObsFilter_3D', 'data'),
  Input('STORE_cellsCtpFilter_3D', 'data'),
  Input('STORE_ifmulti_3D', 'data'),
  Input('BUTTON_singlePlot_3D', 'n_clicks'),
  Input('BUTTON_multiPlot_3D', 'n_clicks'),
  State('DROPDOWN_singleName_3D', 'value'),
  State('STORE_multiNameInfo_3D', 'data'),
  backgroud = True,
  manager = background_callback_manager
)
def update_spatial_plotFeature3D_ctpExpViolin(featureType, stage, cellsExp, cellsObs, cellsCtp, ifmulti, splot, mplot, sname, minfo):
  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]
  
  cells_to_use = list(set(cellsObs) & set(cellsExp) & set(cellsCtp))
  adata = adata[cells_to_use]

  if not ifmulti:
    fig = show_ctpExpViolin(adata, sname)
  else:
    fig = show_multiFeatures_ctpExpViolin(adata, minfo)
    
  return fig

# moran SVG offcanvas
@app.callback(
  Output('OFFCANVAS_moranRes_3D', 'is_open'),
  Input('BUTTON_showMoran_3D', 'n_clicks'),
  prevent_initial_call = True
)
def show_moranRes_offcanvas(click):
  if click:
    return True
  
@app.callback(
  Output('DATATABLE_moranRes_3D', 'data'),
  Output('DATATABLE_moranRes_3D', 'columns'),
  
  Input('BUTTON_calMoran_3D', 'n_clicks'),
  State('STORE_cellsCtpFilter_3D', 'data'),
  State('STORE_cellsExpFilter_3D', 'data'),
  State('STORE_cellsObsFilter_3D', 'data'),
  State('DROPDOWN_stage_3D', 'value'),
  State('DROPDOWN_featureType_3D', 'value'),
  prevent_initial_call=True,
  background = True,
  manager = background_callback_manager,
  running = [
    (Output('BUTTON_showMoran_3D', 'disabled'), True, False),
    (Output('BUTTON_calMoran_3D', 'children'), 'Waiting...\n(<1min)', 'Calculate\n(SVGs)'),
    (Output('BUTTON_calMoran_3D', 'color'), 'danger', 'dark'),
    (Output('BUTTON_calMoran_3D', 'disabled'), True, False),
    (Output('OFFCANVAS_moranRes_3D', 'is_open'), False, True),
  ]
)
def cal_moranRes(click, filter1, filter2, filter3, stage, featureType):
  
  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_data[stage]

  cells = list( set(filter1) & set(filter2) & set(filter3) )
  
  df = cal_moran_3D(adata[cells])
  df = df.reset_index(names='Feature')
  return (df.to_dict('records'),
          [
            {"name": i, "id": i, "deletable": False, 'type': 'numeric', 
              'format':Format(precision=4)} 
            for i in df.columns
          ]
        )

# In[] run:

tabs = dbc.Col(
  spatial_tabs,
  id = 'tabs'
)

app.layout = dbc.Container(
  [
    header,
    dbc.Row([
      dbc.Col([
        tabs,
      ], width=12)
    ],)
  ],
  fluid=True,
)

if __name__ == '__main__':
  app.run_server(
    host='10.193.0.208',
    port='8052',
    debug=True,
    jupyter_mode = 'external'
  )
  

# In[] test set/bool speed

adata = exp_data['E7.75']

cells = adata.obs_names.to_list()
a = np.random.choice(cells, size=22000, replace=False)
b = np.random.choice(cells, size=22000, replace=False)
c = np.random.choice(cells, size=22000, replace=False)
d = np.random.choice([True, False], size=len(cells))
e = np.random.choice([True, False], size=len(cells))
f = np.random.choice([True, False], size=len(cells))


# %%
adata = exp_data['E7.75']
tmp = list(set(a) & set(b) & set(c))
adata = adata[tmp]
# %%
adata = exp_data['E7.75']
adata = adata[d&e&f]

# %%
stage = 'E7.75'
featureType = 'Gene'
cellsObs = a
cellsExp = b
cellsCtp = c
ifmulti = False
sname = 'T'

if featureType == 'Gene':
    adata = exp_data[stage]
elif featureType == 'Regulon':
    adata = auc_data[stage]

cells_to_use = list(set(cellsObs) & set(cellsExp) & set(cellsCtp))
adata = adata[cells_to_use]

if not ifmulti:
  fig = show_expViolin(adata, sname)