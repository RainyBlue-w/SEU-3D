#!/usr/bin/env python
# coding: utf-8

# # prepare

# ## env

# In[254]:


import dash
dash.register_page(__name__)

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
import loompy as lp
import math

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt
import re
import seaborn as sns
from concurrent import futures

import diskcache



background_callback_manager = DiskcacheManager(diskcache.Cache("/rad/wuc/dash_data/reik/cache"))


# ## data

# In[2]:


colors = pd.read_csv("/rad/wuc/Reik/celltype_color.csv")
colors['celltype'] = [re.sub(' |/', '_', i) for i in colors['celltype']]
celltype_colors = dict(zip(colors['celltype'], colors['colour']))


# Code below shows how the anndata with diffgene results was generated:

# In[23]:


# prepare data
stage_list = ['E7.5','E7.75','E8.0','E8.5','E8.75']

umap = pd.read_csv('/rad/wuc/Reik/all_combine_umap.csv', index_col=0)
exp_data = sc.read_h5ad('/rad/wuc/Reik/anndata.h5ad')
exp_data = exp_data[(exp_data.obs['sample'] != 'E8.5_CRISPR_T_WT')&(exp_data.obs['sample'] != 'E8.5_CRISPR_T_KO')]
sc.pp.normalize_total(exp_data, target_sum=1e4)
sc.pp.log1p(exp_data)
all(exp_data.obs.index.isin(umap.index))


# In[25]:


tmp = {}
for stage in stage_list:
  tmp[stage] = exp_data[exp_data.obs['stage'] == stage]
exp_data = tmp.copy()
del(tmp)

tmp = {}
for stage in stage_list:
  tmp[stage] = umap[umap.index.str.startswith(stage)]
umap = tmp
del(tmp)


# In[10]:


auc_mtx = {}
regulon_geneset = {}
for stage in stage_list:
  fpath = f"/rad/wuc/Reik/scenic/Reik_{stage}_pyscenic_output.loom"
  lf = lp.connect( fpath, mode='r+', validate=False )
  auc_mtx[stage] = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
  df = pd.DataFrame(lf.ra.Regulons, index = lf.ra.Gene).astype("bool")
  regulon_geneset[stage] = {}
  for regulon in df.columns:
    regulon_geneset[stage][regulon] = lf.ra.Gene[df.loc[:,regulon]].tolist()
  lf.close()

auc_data = {}
for stage in stage_list:
  auc_data[stage] = sc.AnnData(X=auc_mtx[stage], obs=exp_data[stage].obs)
del(auc_mtx)


# ## plot functions

# In[378]:


def show_feature_umap(adata, feature, embedding, cmap = None, sort=False, ascending=True, **kws):
    if cmap is None:
        cmap = [(0.00, "rgb(244,244,244)"),
                (0.05, "rgb(244, 244, 244)"),
                (1.00, "rgb(34, 94, 168)")
                ]
    embedding.columns = ['umap_1', 'umap_2']
    pdf = pd.concat([embedding, adata[:,feature].to_df()], axis=1)
    if sort is True:
      pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter(
    		data_frame = pdf,
    		x = 'umap_1', y = 'umap_2', color = feature,
    		color_continuous_scale = cmap,
        **kws
    	)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=4,
                      marker_opacity=1)
    plot.update_layout(
      margin=dict(l=0, r=0, t=0, b=0),
      plot_bgcolor = '#ffffff', 
      uirevision='constant',
      coloraxis = {
        'colorbar' : {
          'tickformat': '4.2f',
          'orientation': 'h',
          'y': -0.1
        }
      },
      legend_itemsizing = 'constant',
    )
    plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                               font_size = 20)) 
    return plot


# In[379]:


def show_celltype_umap(adata, embedding, cmap = celltype_colors, sort=False, **kws):
  embedding.columns = ['umap_1', 'umap_2']
  pdf = pd.concat([embedding, adata.obs[['celltype']]], axis=1)
  if sort:
    pdf = pdf.sort_values(by='celltype')
  plot = px.scatter(
  	data_frame = pdf,
    x = 'umap_1', y = 'umap_2', color = 'celltype',
    color_discrete_map = cmap,
    **kws
  )
  plot.update_yaxes(visible=False)
  plot.update_xaxes(visible=False)
  plot.update_traces(marker_size=3,
                    marker_opacity=1)
  plot.update_layout(
    margin=dict(l=0, r=0, t=0, b=0),
    plot_bgcolor = '#ffffff',
    title = '',
    legend = dict(
        title = '',
        orientation='h',
        yanchor='top',
        xanchor='center',
        y=-0.05,
        x=0.5,
    ),
    uirevision='constant',
    legend_itemsizing = 'constant'
  )
  plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                      font_size = 20)) 

  return plot


# # app

# In[380]:




# ## page-Reik

# ## widgets

# In[382]:


reik_dropdown_plotFeature_gene = html.Div(
  [
    dbc.Label("Gene to display"),
    dcc.Dropdown(
      exp_data['E7.5'].var_names,
      'Hand1',
      id="reik_dropdown_plotFeature_gene",
      clearable=True,
      searchable=True,
    ),
  ],
  className="mb-4",
)

reik_dropdown_plotFeature_regulon = html.Div(
  [
    dbc.Label("Regulon to display"),
    dcc.Dropdown(
      auc_data['E7.5'].var_names,
      'Hand1(+)',
      id="reik_dropdown_plotFeature_regulon",
      clearable=True,
      searchable=True,
    ),
  ],
  className="mb-4",
)

reik_table_regulonTargetGene = html.Div(
  [
    dbc.Label("Regulon's target gene"),
    dash_table.DataTable(sort_action="native", page_action='native',
                         page_current= 0, page_size= 10,
                         id = "reik_plotFeature_regulonTargetGene",fill_width=True,
                         style_table={'overflowX': 'auto'},
                         style_cell={
                           'padding-right': '30px',
                           'padding-left': '10px',
                           'text-align': 'center',
                           'marginLeft': 'auto',
                           'marginRight': 'auto'
                         },
                         filter_action="native", 
                        ),
  ],
  className="mb-4",
)


# In[383]:


reik_celltypeGraph = dcc.Graph(figure={}, id="reik_plotFeature_celltypeGraph",style={'height': "75vh"})

reik_featureGraphs = dbc.Row([
  dbc.Col([
    html.H5('Gene', style={'textAlign':'center'}, id='reik_plotFeature_geneTitle'),
    dcc.Graph(figure={}, id='reik_plotFeature_geneGraph', style={'height': '60vh'})
  ], width=12, xxl=6),
  dbc.Col([
    html.H5('Regulon', style={'textAlign':'center'}, id='reik_plotFeature_regulonTitle'),
    dcc.Graph(figure={}, id='reik_plotFeature_regulonGraph', style={'height': '60vh'})
  ], width=12, xxl=6),
])


# ### tabs

# In[384]:


num_stages = [ float(re.sub('E','',i)) for i in stage_list]
num_stages = [int(i) if (i % 1 == 0) else i for i in num_stages ]


# In[385]:


reik_tab_plotFeature = dbc.Tab(
  [dbc.Row([
    # control panel
    dbc.Col([
      dbc.Card(
        [
          reik_dropdown_plotFeature_gene,
          reik_dropdown_plotFeature_regulon,
          reik_table_regulonTargetGene,
        ],
        body=True,
        id = 'reik_plotFeature_controlPanel'
      )
    ], width=2),
    
    # fig panel
    dbc.Col([
      dbc.Row([
        dbc.Col([
          html.H5('celltype', style={'textAlign':'center'}),
          reik_celltypeGraph,
        ], width=4),
        dbc.Col([
          reik_featureGraphs,
        ], width=8),
        dbc.Col([
          dcc.Slider(id='reik_plotFeature_stageSlider',
            min = math.floor(num_stages[0]),
            max = math.ceil(num_stages[-1]),
            value=num_stages[0], step=None,
            marks= dict(zip(
              num_stages,
              stage_list))
          ),
        ], width=12),
      ],)
    ], width=10)
    
  ])],
  label = "Plot feature",
  tab_id = "reik_tab_plotFeature"
)


# In[386]:


reik_tabs = dbc.Card(
    dbc.Tabs(
        [reik_tab_plotFeature],
        active_tab = "reik_tab_plotFeature",  
        id = "reik_tabs",
    ),
)


# ### callbacks

# In[387]:


@callback(
  Output('reik_dropdown_plotFeature_gene', 'options'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_gene', 'search_value')
)
def update_dropdown_options_gene(stage, search):
  if not stage:
    raise PreventUpdate
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  vars = exp_data[stage].var_names
  return vars[vars.str.startswith(search)].sort_values()

@callback(
  Output('reik_dropdown_plotFeature_regulon', 'options'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_regulon', 'search_value')
)
def update_dropdown_options_regulon(stage, search):
  if not stage:
    raise PreventUpdate
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  vars = auc_data[stage].var_names
  return vars[vars.str.startswith(search)].sort_values()


# In[388]:


@callback(
  Output('reik_plotFeature_celltypeGraph','figure'),
  Input('reik_plotFeature_stageSlider', 'value'),
)
def update_celltypeGraph(stage):
  if not stage:
    raise PreventUpdate
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  return show_celltype_umap(adata=exp_data[stage], embedding=umap[stage], cmap=celltype_colors, sort=True)


@callback(
  Output('reik_plotFeature_geneGraph','figure'),
  Output('reik_plotFeature_geneTitle', 'children'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_gene', 'value')
)
def update_geneGraph(stage, gene):
  if not stage or not gene:
    raise PreventUpdate
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  return (
    show_feature_umap(adata=exp_data[stage], feature=gene, embedding=umap[stage], sort=True),
    gene
  )

@callback(
  Output('reik_plotFeature_regulonGraph','figure'),
  Output('reik_plotFeature_regulonTitle', 'children'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_regulon', 'value')
)
def update_regulonGraph(stage, regulon):
  if not stage or not regulon:
    raise PreventUpdate
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  return (
    show_feature_umap(adata=auc_data[stage], feature=regulon, embedding=umap[stage], sort=True),
    regulon
  )

@callback(
  Output('reik_plotFeature_regulonTargetGene', 'data'),
  Output('reik_plotFeature_regulonTargetGene', 'columns'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_regulon', 'value')
)
def update_regulonTargetGenes(stage, regulon):
  if not stage or not regulon:
    raise PreventUpdate
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  df = pd.DataFrame(regulon_geneset[stage][regulon],
            columns=[f"{regulon}'s target gene"])
  df = df.reset_index().rename(columns={"index": "id"})
  return df.to_dict('records'), [{"name": i, "id": i, "deletable": False} for i in df.columns if i != 'id']

@callback(
  Output('reik_dropdown_plotFeature_gene', 'value'),
  Input('reik_plotFeature_regulonTargetGene', 'active_cell'),
  Input('reik_plotFeature_stageSlider', 'value'),
  State('reik_dropdown_plotFeature_regulon', 'value')
)
def update_regulonGraph_by_targetGene(active_cell, stage, regulon):
  if not active_cell or not regulon:
    raise PreventUpdate
  if not stage:
    raise PreventUpdate
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  row = active_cell['row_id']
  df = pd.DataFrame(regulon_geneset[stage][regulon],
        columns=[f"{regulon}'s target gene"])
  return df.iloc[row,0]


# ## all layout

# In[389]:


tabs = dbc.Col(
  reik_tabs,
  id = 'tabs'
)
layout = dbc.Container(
  [
    dbc.Row([
      dbc.Col([
        tabs,
      ], width=12)
    ],)
  ],
  fluid=True,
  className="Container-all",
)


# ## run

# In[391]:



# In[ ]:





# In[ ]:




