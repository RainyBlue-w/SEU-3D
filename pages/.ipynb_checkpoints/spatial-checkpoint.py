#!/usr/bin/env python
# coding: utf-8

# # env

# In[1]:


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


# In[2]:


os.getpid()


# # data

# In[3]:


exp_data = {
  'E7.5': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.5_HC0.5_min400.h5ad"),
  'E7.75': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.75_HC0.5_min400.h5ad"),
  'E8.0': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E8.0_HC0.5_min400.h5ad")
}


# In[4]:


for stage in exp_data.keys():
  exp_data[stage].raw = exp_data[stage].copy()
  sc.pp.normalize_total(exp_data[stage], target_sum=1e4)
  sc.pp.log1p(exp_data[stage])


# In[5]:


for k,v in exp_data.items():
  v.obs.germ_layer = [i.replace('exe-ecto', 'ectoderm') for i in v.obs.germ_layer.values]


# In[6]:


auc_data = {}
regulon_geneset = {}
for stage in ['E7.5', 'E7.75', 'E8.0']:
  h5 = h5py.File( '/rad/wuc/dash_data/spatial/matrix_data/%s_auc_mtx.h5' % (stage))
  auc_mtx = pd.DataFrame(h5['auc_matrix'], index=h5['cell'][:].astype(str), columns=h5['regulon'][:].astype(str))
  auc_data[stage] = sc.AnnData(X=auc_mtx, obs=exp_data[stage].obs)
  regulon_geneset[stage] = json.loads(h5['geneset'][()])
  h5.close()
del(auc_mtx)


# In[7]:


genes_min_pval = pd.read_csv("/rad/wuc/dash_data/spatial/sparkX_res/genes_padj_minVal.csv",
                             header=[0], index_col=[0,1]) 


# In[8]:


genes_all_pval = pd.read_csv("/rad/wuc/dash_data/spatial/sparkX_res/genes_padj_combine.csv",
                             header=[0,1], index_col=[0,1]) 
genes_all_pval = genes_all_pval.loc[:,[('ecto', 'adjustedPval'), ('meso', 'adjustedPval'), ('endo', 'adjustedPval'), ('all', 'adjustedPval')]]
genes_all_pval.columns = ['ecto p.adj', 'meso p.adj', 'endo p.adj', 'all p.adj']
genes_all_pval = genes_all_pval.groupby(level=0, group_keys=False
                                               ).apply(lambda x: x.sort_values(by='all p.adj'))


# In[9]:


ctp_cmap = pd.read_csv("/rad/wuc/dash_data/spatial/celltype_cmap.csv")
ctp_cmap = dict(zip(ctp_cmap['celltype'], ctp_cmap['Epiblast']))


# In[10]:


# data_for_dask = {}
# data_for_dask['exp_data'] = {}
# for stage in exp_data.keys():
#   data_for_dask['exp_data'][stage] = 


# # plot functions

# In[11]:


def show_celltype_spatial(adata, embedding, cmap = ctp_cmap, **kws):
  pdf = pd.DataFrame(np.array(embedding), 
                    index=adata.obs_names, 
                    columns=['ebd1', 'ebd2'])
  pdf = pd.concat([pdf, adata.obs[['celltype','germ_layer']]], axis=1)
  pdf = pdf.sort_values(by='celltype')
  plot = px.scatter(
  	data_frame = pdf,
    x = 'ebd1', y = 'ebd2', color = 'celltype',
    color_discrete_map = cmap,
    **kws
  )
  plot.update_yaxes(visible=False)
  plot.update_xaxes(visible=False)
  plot.update_traces(marker_size=4.5,
                    marker_opacity=1)
  plot.update_layout(
    margin=dict(l=10, r=10, t=30, b=0),
    plot_bgcolor = '#ffffff', 
    uirevision='constant',
    legend_itemsizing = 'constant'
  )
  plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                      font_size = 20)) 

  return plot


# In[12]:


def show_feature_spatial(adata, feature, embedding, cmap = None, sort=False, ascending=True, **kws):
    if cmap is None:
        cmap = [(0.00, "rgb(244,244,244)"),
                (0.05, "rgb(244, 244, 244)"),
                (1.00, "rgb(34, 94, 168)")
                ]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['ebd1', 'ebd2'])
    pdf = pd.concat([pdf, adata[:,feature].to_df(), adata.obs], axis=1)
    if sort is True:
      pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter(
    		data_frame = pdf,
    		x = 'ebd1', y = 'ebd2', color = feature,
    		color_continuous_scale = cmap,
        **kws
    	)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=4.5,
                      marker_opacity=1)
    plot.update_layout(
      margin=dict(l=10, r=10, t=30, b=0),
      plot_bgcolor = '#ffffff', 
      uirevision='constant',
      coloraxis = {
        'colorbar' : {'tickformat': '4.2f'}
      }
    )
    plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                               font_size = 20)) 
    return plot


# In[13]:


def show_feature_spatial(adata, feature, embedding, cmap = None, sort=False, ascending=True, **kws):
    if cmap is None:
        cmap = [(0.00, "rgb(244,244,244)"),
                (0.05, "rgb(244, 244, 244)"),
                (1.00, "rgb(34, 94, 168)")
                ]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['ebd1', 'ebd2'])
    pdf = pd.concat([pdf, adata[:,feature].to_df(), adata.obs], axis=1)
    if sort is True:
      pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter(
    		data_frame = pdf,
    		x = 'ebd1', y = 'ebd2', color = feature,
    		color_continuous_scale = cmap,
        **kws
    	)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=4.5,
                      marker_opacity=1)
    plot.update_layout(
      margin=dict(l=10, r=10, t=30, b=0),
      plot_bgcolor = '#ffffff', 
      uirevision='constant',
      coloraxis = {
        'colorbar' : {'tickformat': '4.2f'}
      }
    )
    plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                               font_size = 20)) 
    return plot


# In[14]:


def show_celltype_spatial_3D(adata, embedding, cmap = ctp_cmap, **kws):

    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['ebd1', 'ebd2', 'ebd3'])
    pdf = pd.concat([pdf, adata.obs['celltype']], axis=1)
    plot = px.scatter_3d(
    		data_frame = pdf,
    		x = 'ebd1', y = 'ebd2', z='ebd3',
        color = 'celltype',
    		color_discrete_map=cmap,
        **kws
    	)
    plot.update_traces(marker_size=3,
                      marker_opacity = 1)
    plot.update_layout(
      margin=dict(l=10, r=10, t=30, b=0),
      uirevision='constant',
      legend_itemsizing = 'constant',
      scene = dict(
        xaxis = {'title': 'x',
              'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        yaxis = {'title': 'y',
              'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        zaxis = {'title': 'z',
              'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        bgcolor = 'white'
      )
    )
    return plot


# In[15]:


def show_feature_spatial_3D(adata, feature, embedding, cmap = None, sort=False, ascending=True, **kws):
    if cmap is None:
        cmap = [(0.00, "rgb(244,244,244)"),
                (0.05, "rgb(244, 244, 244)"),
                (1.00, "rgb(34, 94, 168)")
                ]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['ebd1', 'ebd2', 'ebd3'])
    pdf = pd.concat([pdf, adata[:,feature].to_df()], axis=1)
    if sort is True:
      pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter_3d(
    		data_frame = pdf,
    		x = 'ebd1', y = 'ebd2', z='ebd3',
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
      },
      scene = dict(
        xaxis = {'title': 'x',
              'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        yaxis = {'title': 'y',
              'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        zaxis = {'title': 'z',
              'visible': True,
                'backgroundcolor' :"white",
                'showbackground': True,
                'zerolinecolor': 'grey',
                'gridcolor': 'grey',
                'nticks': 6},
        bgcolor = 'white'
      )
    )

    plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                               font_size = 20))
    return plot


# In[16]:


def show_features_spatial_regularExp(adata, stage, pattern, odir, featureType, embedding, cmap=None, sort=False, ascending=True, dpi=100, **kws):
  
    img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features = [i for i in features if i in genes_min_pval.loc[stage].index.to_list()]
    ordered_features = genes_min_pval.loc[stage].loc[features].sort_values(by='PadjMinVal', ascending=True).index.to_list()
  
    if os.path.exists( img_dir ):
      return img_dir
      
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['ebd1', 'ebd2'])
    
    features_df = adata[:,ordered_features].to_df()
    features_df.columns = ['%s\n(%.2e)' % (i,genes_min_pval.loc[stage].loc[i, 'PadjMinVal']) for i in features_df.columns]

    pdf = pd.concat([pdf, features_df, adata.obs.germ_layer], axis=1)
    pdf = pd.melt(pdf, id_vars = ['ebd1', 'ebd2', 'germ_layer'], var_name='feature', value_name='value')
    if sort:
      pdf = pdf.sort_values(by='value', ascending=ascending)
    pdf['feature'] = pdf['feature'].astype('category').values.reorder_categories(features_df.columns)
    pdf['germ_layer'] = pdf['germ_layer'].astype('category').values.reorder_categories(['ectoderm','mesoderm','endoderm'])

    # plotnine
    (
      ggplot() + 
      geom_point(data = pdf, mapping=aes(x='ebd1', y='ebd2', color='value')) +
      facet_grid('feature ~ germ_layer') + 
      scale_color_gradientn(colors = ['#e0e0e0', '#e0e0e0','#225ea8'], values = [0,0.05,1]) +
      theme(
        # legend_position='top',
        # legend_direction='horizontal',
        axis_title = element_text(size=16),
        axis_text = element_text(size=10),
        panel_grid = element_blank(),
        panel_border = element_rect(linewidth=0.4, fill= 'None'),
        panel_background = element_rect(fill= 'None'),
        strip_text_x = element_text(size=16),
        strip_text_y = element_text(size=16, face = 'bold', angle=0)
      )
    ).save(img_dir, width=12, height=1+3*len(features), dpi=dpi, 
           limitsize=False, verbose=False)
    return img_dir


# In[17]:


def show_featuresCtpcounts_spatial_regularExp(adata, stage, pattern, odir, featureType, embedding, cmap=ctp_cmap, sort=False, ascending=True, dpi=150, **kws):
  
  img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/plotFeatureSeries_%s_ctpCounts_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
  
  if os.path.exists( img_dir ):
    return img_dir
  
  features = [i  for i in adata.var_names if re.match(pattern, i)]
  features = [i for i in features if i in genes_min_pval.loc[stage].index.to_list()]
  ordered_features = genes_min_pval.loc[stage].loc[features].sort_values(by='PadjMinVal', ascending=True).index.to_list()
  
  ctp_counts = {}
  for gene in ordered_features:
    df = adata[:,gene].to_df()
    # thr = df.min()+(df.max()-df.min())*0.05
    df = df[df[gene] > 0]
    counts = pd.DataFrame(adata.obs['celltype'][df.index].value_counts())
    counts['gene'] = gene
    counts['count'] = counts['count']/sum(counts['count'])
    ctp_counts[gene] = counts
  ctp_counts = pd.concat(ctp_counts, axis=0)
  ctp_counts['celltype'] = np.array(ctp_counts.index.to_list())[:,1]
  ctp_counts['text_y'] = ctp_counts['count'].max()/2
  ctp_counts['gene'] = ctp_counts['gene'].astype('category').values.reorder_categories(ordered_features)
  bar_df = ctp_counts[ctp_counts['count']>0].groupby(['gene']).head(5)
  bar_df['ctp_order'] = bar_df['celltype'] +'_' + bar_df['gene'].astype(str)
  bar_df['ctp_order'] = bar_df['ctp_order'].astype('category').values.reorder_categories(bar_df['ctp_order'][::-1])
  (
    ggplot() + 
      geom_bar(data = bar_df, mapping=aes(x='ctp_order', y='count', fill = 'celltype'), stat='summary') +
      facet_grid(
        'gene~.', scales='free_y'
                ) +
      scale_y_continuous( labels = lambda list: [("%.2f%%" % (x*100)) for x in list] ) +
      scale_fill_manual(breaks=list(cmap.keys()),
                        values=list(cmap.values())) +
      guides(
        fill= guide_legend(
          ncol = 1,
          byrow = True
        )
      ) + 
      theme(legend_position='none',
        axis_text_y = element_text(size=12),
        axis_text_x = element_text(size=12),
        axis_title_x = element_blank(),
        axis_ticks_major_x = element_blank(),
        axis_ticks_minor_x = element_blank(),
        axis_title_y=element_blank(),
        strip_text_y = element_text(size=16, face = 'bold', angle=0)
      ) +
      scale_x_discrete(labels = lambda list: [re.search(r'^(.*?)_',x).group(1) for x in list]) + 
      # geom_text(data = ctp_counts[ctp_counts['count'] > 0.05],
      #           mapping = aes(x='celltype', y='text_y', label = 'celltype'),
      #           size=8,
      #           angle=90) + 
      coord_flip()
  ).save(img_dir, width=8, height=1+3*len(ordered_features), dpi=dpi, 
           limitsize=False, verbose=False)
  return img_dir


# In[18]:


def show_metadata_discrete_spatial_germLayer(adata, stage,  odir, embedding, obs=['celltype'], cmap = ctp_cmap, dpi=150, sort = False, ascending=True, **kws):
  
  img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/poltFeatureSeries_%s_plot_dpi%d.png' % (odir, stage, obs, dpi)

  # if os.path.exists( img_dir ):
  #   return img_dir

  pdf = pd.DataFrame(np.array(embedding), 
                     index=adata.obs_names, 
                     columns=['ebd1', 'ebd2'])
  
  obs_df = pd.DataFrame(adata.obs[obs])

  pdf = pd.concat([pdf, obs_df, adata.obs.germ_layer], axis=1)
  pdf = pd.melt(pdf, id_vars = ['ebd1', 'ebd2', 'germ_layer'], var_name='obs', value_name='value')
  if sort:
    pdf = pdf.sort_values(by='value', ascending=ascending)
  pdf['obs'] = pdf['obs'].astype('category').values.reorder_categories(obs_df.columns)

  # plotnine
  (
    ggplot() + 
    geom_point(data = pdf, mapping=aes(x='ebd1', y='ebd2', color='value'), size=0.4) +
    facet_grid('obs ~ germ_layer') + 
    # scale_color_gradientn(colors = ['#e0e0e0', '#e0e0e0','#225ea8'], values = [0,0.05,1]) +
    scale_color_manual(values = ctp_cmap) +
    theme(
      # legend_position='top',
      # legend_direction='horizontal',
      axis_title = element_text(size=16),
      axis_text = element_text(size=10),
      panel_grid = element_blank(),
      panel_border = element_rect(linewidth=0.4, fill= 'None'),
      panel_background = element_rect(fill= 'None'),
      strip_text_x = element_text(size=16),
      strip_text_y = element_text(size=16, face = 'bold', angle=0)
    )
  ).save(img_dir, width=12, height=1+3*len(obs), dpi=dpi, 
         limitsize=False, verbose=False)
  return img_dir


# In[19]:


fig1 =  show_feature_spatial(exp_data['E7.75'], 'Cdx1', embedding = exp_data['E7.75'].obs[['x_flatten', 'y_flatten']], sort=True, ascending=True,
                            facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']})
fig2 =  show_feature_spatial_3D(exp_data['E7.75'], 'Cdx1', embedding = exp_data['E7.75'].obs[['x', 'y', 'z']], sort=True, ascending=True)


# # app layout

# In[288]:


dbc_css = "/home/wuc/dashapps/css/dbc.min.css"
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc_css])


# ## header

# In[289]:


header = dbc.NavbarSimple(
    [
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("dataset", header=True),
                dbc.DropdownMenuItem("atlas"),
                dbc.DropdownMenuItem("spatial"),
            ],
            nav=True,
            in_navbar=True,
            label="More",
        ),
    ],
    brand="Omics-viewer",
    color="dark",
    dark=True,
    sticky='top',
    style = {"height": "6vh"}
)


# ## page-atlas

# ## page-spatial

# ### dropdowns

# In[290]:


spatial_dropdown_featureType = html.Div(
    [
        dbc.Label("Feature type"),
        dcc.Dropdown(
            ['Gene', 'Regulon'],
            'Gene',
            id="spatial_dropdown_featureType",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)
spatial_dropdown_featureName = html.Div(
    [
        dbc.Label("Feature name"),
        dcc.Dropdown(
            exp_data['E7.75'].var_names,
            'Cdx1',
            id="spatial_dropdown_featureName",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

spatial_dropdown_plotDim = html.Div(
    [
        dbc.Label("Plot type"),
        dcc.Dropdown(
            ['2D', '3D'],
            '2D',
            id="spatial_dropdown_plotDim",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

spatial_dropdown_featureType_series = html.Div(
    [
        dbc.Label("Feature type"),
        dcc.Dropdown(
            ['Gene'],
            'Gene',
            id="spatial_dropdown_featureType_series",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)
spatial_input_featureName_series = html.Div(
    [
        dbc.Label("Starts with"),
        dbc.InputGroup(
          [
            dbc.Input(id="spatial_input_featureName_series"),
            dbc.Button('Plot', id='spatial_inputButton_featureName_series_plot',
                                n_clicks=0, color='primary'),
          ]
        )
    ],
    className="mb-4",
)
spatial_dropdown_stage = html.Div(
    [
        dbc.Label("Stage"),
        dcc.Dropdown(
            ['E7.5', 'E7.75', 'E8.0'],
            'E7.75',
            id="spatial_dropdown_stage",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

spatial_dropdown_stage_series = html.Div(
    [
        dbc.Label("Stage"),
        dcc.Dropdown(
            ['E7.5', 'E7.75', 'E8.0'],
            'E7.75',
            id="spatial_dropdown_stage_series",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)


# ### tabs

# In[291]:


spatial_tab_plotFeature = dbc.Tab(
  [dbc.Row([
    dbc.Col([
      dbc.Card(
        [
          spatial_dropdown_featureType,
          spatial_dropdown_featureName,
          spatial_dropdown_plotDim,
          spatial_dropdown_stage
        ],
        body=True,
        id = 'spatial_control_plotFeature'
      ),      
    ], width=2),
    dbc.Col([dbc.Row([
      dbc.Col([
        dcc.Graph(figure={}, id="spatial_plotFeature_graph", style={'height': "45vh"}),
      ],align = "center",className="g-0", width=12,id='spatial_plotFeature_graph_col'),
      dbc.Col([
        dcc.Graph(figure={}, id="spatial_plotFeature_graph_ctp", style={'height': "45vh"}),
      ],align = "center",className="g-0", width=12,id='spatial_plotFeature_graph_ctp_col')
    ]),],width=10),
  ])],
  label = "Plot feature",
  tab_id = "spatial_tab_plotFeature"
)

spatial_tab_plotFeatureSeries = dbc.Tab(
  [dbc.Row([
    dbc.Col([
      dbc.Card(
        [
            spatial_dropdown_featureType_series,
            spatial_input_featureName_series,
            spatial_dropdown_stage_series,
        ],
        body=True,
        id = 'spatial_control_plotFeatureSeries'
      )
    ], width=2),
    dbc.Col([
      dbc.Row(
        [
          dbc.Col([
            dcc.Graph(figure={}, id="spatial_plotFeatureSeries_graph_ctp", style={'height': "50vh"}),
          ],align = "center",className="g-0", width=8),
        ]
      ),
      dbc.Row([
        dbc.Col(
          [
            html.Img(
              src = Image.open('/rad/wuc/dash_data/spatial/tmp/plotFeatureSeries/E7.75/plotFeatureSeries_E7.75_^Cdx_Gene_dpi100.png'),
                     id = 'spatial_plotFeatureSeries_img', style = {'width': '90vh'})
          ],
          align = "center", className="g-0", width=7),
        dbc.Col(
          [
            html.Img(
              src = Image.open('/rad/wuc/dash_data/spatial/tmp/plotFeatureSeries/E7.75/plotFeatureSeries_E7.75_ctpCounts_^Cdx_Gene_dpi150.png'),
                     id = 'spatial_plotFeatureSeries_ctpCounts_img', style = {'width': '60vh'})
          ],
          
          align = "center", className="g-0", width=5)
      ],),
    ], width=10)
  ])],
  label = "Plot feature(Series)",
  tab_id = "spatial_tab_plotFeatureSeries"
)

spatial_tab_plotFeatureSparkx = dbc.Tab(
  [
    dcc.Store(
      id='clientside_store_table_sparkx'
    ),
    dbc.Row([
      dbc.Col([
        dbc.Card(
          dbc.Row([
            dbc.Col([
              dbc.Label("Feature type"),
              dcc.Dropdown(['Gene'], 'Gene', id='spatial_dropdown_featureType_sparkx', clearable=False)
            ], width=4),
            dbc.Col([
              dbc.Label("Stage"),
              dcc.Dropdown(['E7.5', 'E7.75', 'E8.0'],'E7.75', id='spatial_dropdown_stage_sparkx', clearable=False)
            ], width=4),
            dbc.Col([
              dbc.Label("Genes per page"),
              dcc.Dropdown([10, 15, 20, 30, 50, 100], 15, id='spatial_dropdown_pageSize_sparkx', clearable=False)
            ], width=4)
          ])
        ),
        dash_table.DataTable(id = "spatial_dataTable_sparkx",
          row_selectable = 'single', sort_action="native", page_action='native',
          filter_action="native", page_current= 0, page_size= 15,fill_width=True, 
          style_table={'overflowX': 'auto'},
          style_cell={'padding-right': '10px', 'padding-left': '10px',
            'text-align': 'center', 'marginLeft': 'auto', 'marginRight': 'auto' }
        ),
      ], width=4),
      dbc.Col([
        dcc.Graph(figure=fig1, id="spatial_featureGraph_sparkx", style={'height': "36vh"}),
        dcc.Graph(figure={}, id="spatial_celltypeGraph_sparkx", style={'height': "36vh"}),
      ], width=8),
    ]),
  ],
  label = "SPARK-X",
  tab_id = "spatial_tab_plotFeatureSparkx"
)


# ### layout

# In[292]:


# spatial_controls = dbc.Col(
#     spatial_control_plotFeature,
#     id = "spatial_controls",
#     style = {'position': 'fixed'},
#     width=2
# )


# In[293]:


spatial_tabs = dbc.Card(
    dbc.Tabs(
        [spatial_tab_plotFeature, spatial_tab_plotFeatureSeries, spatial_tab_plotFeatureSparkx],
        active_tab = "spatial_tab_plotFeature",  
        id = "spatial_tabs",
    ),
)


# ### callbacks

# In[294]:


@app.callback(
  Output('spatial_dropdown_featureName', 'options'),
  Input('spatial_dropdown_featureName', 'search_value'),
  Input('spatial_dropdown_featureType', 'value'),
  Input('spatial_dropdown_stage', 'value')
)
def update_spatial_dropdown_featureName(search, featureType, stage):
  if not search:
    raise PreventUpdate
  
  if featureType == 'Gene':
    if not search:
      return exp_data[stage].var_names
    else:
      return exp_data[stage].var_names[exp_data[stage].var_names.str.startswith(search)].sort_values()
  elif featureType == 'Regulon':
    if not search:
      return auc_data[stage].var_names
    else:
      return auc_data[stage].var_names[auc_data[stage].var_names.str.startswith(search)].sort_values()

@app.callback(
  Output('spatial_plotFeature_graph', 'figure'),
  Output('spatial_plotFeature_graph', 'style'),
  Output('spatial_plotFeature_graph_col', 'width'),
  State('spatial_dropdown_featureType', 'value'),
  Input('spatial_dropdown_featureName', 'value'),
  Input('spatial_dropdown_plotDim', 'value'),
  Input('spatial_dropdown_stage', 'value'),
  background=True,
  manager=background_callback_manager,
)
def update_spatial_plotFeature_graph(featureType, name, type, stage):
  if name is None:
    raise PreventUpdate
  
  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]
    
  if type == '2D':
    return (
      show_feature_spatial(adata, name, embedding = adata.obs[['x_flatten', 'y_flatten']], sort=True, ascending=True,
                            facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']}),
      {'height': '40vh'}, 12
    )
  elif type == '3D':
    return (
      show_feature_spatial_3D(adata, name, embedding = adata.obs[['x', 'y', 'z']], sort=True, ascending=True),
      {'height': '85vh'}, 5
    )
  else:
    raise PreventUpdate

@app.callback(
  Output('spatial_plotFeature_graph_ctp', 'figure'),
  Output('spatial_plotFeature_graph_ctp', 'style'),
  Output('spatial_plotFeature_graph_ctp_col', 'width'),
  State('spatial_dropdown_featureType', 'value'),
  Input('spatial_dropdown_plotDim', 'value'),
  Input('spatial_dropdown_stage', 'value'),
  background=True,
  manager=background_callback_manager,
)
def update_spatial_plotFeature_ctpGraph(featureType, dim, stage):

  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]
  if dim == '2D':
    return (
      show_celltype_spatial(adata, embedding = adata.obs[['x_flatten', 'y_flatten']],
                           facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']}),
      {'height': '40vh'}, 12
    )
  elif dim == '3D':
    return (
      show_celltype_spatial_3D(adata, embedding = adata.obs[['x', 'y', 'z']]),
      {'height': '85vh'}, 7
    )
  else:
    raise PreventUpdate

@app.callback(
  Output('spatial_plotFeatureSeries_img', 'src'),
  Output('spatial_plotFeatureSeries_ctpCounts_img', 'src'),
  State('spatial_dropdown_featureType_series', 'value'),
  State('spatial_input_featureName_series', 'value'),
  Input('spatial_inputButton_featureName_series_plot', 'n_clicks'),
  State('spatial_dropdown_stage_series', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('spatial_inputButton_featureName_series_plot', 'children'), 'Loading', 'Plot'),
    (Output('spatial_tab_plotFeatureSeries', 'label'), 'Loading...', "Plot feature(Series)"),
    (Output('spatial_inputButton_featureName_series_plot', 'disabled'), True, False),
    (Output('spatial_inputButton_featureName_series_plot', 'color'), 'danger', 'primary'),
    (Output('spatial_input_featureName_series', 'disabled'),True, False),
  ],
  prevent_initial_call = True
)
def update_spatial_plotFeature_graphSeries(featureType, pattern, click, stage):
  if pattern is None:
    raise PreventUpdate

  if click:
    pattern = '^'+pattern
    if featureType == 'Gene':
      adata = exp_data[stage]
    elif featureType == 'Regulon':
      adata = auc_mtx[stage]

    # img_dir1, ordered_features = show_features_spatial_regularExp(adata, stage, pattern, 'plotFeatureSeries', featureType, embedding = adata.obs[['x_flatten', 'y_flatten']], sort=True, ascending=True,)
    # img_dir2 = show_featuresCtpcounts_spatial_regularExp(adata, stage, pattern, ordered_features,'plotFeatureSeries', featureType, embedding = adata.obs[['x_flatten', 'y_flatten']], sort=True, ascending=True,)
    with futures.ThreadPoolExecutor(max_workers=8) as executor:
      t1 = executor.submit(show_features_spatial_regularExp, adata, stage, pattern, 'plotFeatureSeries', featureType, embedding = adata.obs[['x_flatten', 'y_flatten']], sort=True, ascending=True)
      t2 = executor.submit(show_featuresCtpcounts_spatial_regularExp, adata, stage, pattern,'plotFeatureSeries', featureType, embedding = adata.obs[['x_flatten', 'y_flatten']], sort=True, ascending=True)
      img_dir1 = t1.result()
      img_dir2 = t2.result()
    return Image.open(img_dir1), Image.open(img_dir2)
  else:
    raise PreventUpdate


# In[295]:


@app.callback(
  Output('spatial_plotFeatureSeries_graph_ctp', 'figure'),
  Input('spatial_dropdown_stage_series', 'value'),
)
def update_spatial_plotFeatureSeries_ctpGraph(stage):
  fig = show_celltype_spatial(exp_data[stage], embedding = exp_data[stage].obs[['x_flatten', 'y_flatten']],
                             facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']})
  fig.update_layout(
    legend = {
      'orientation': 'h'
    }
  )
  return fig


# In[300]:


@app.callback(
  [
    Output('spatial_dataTable_sparkx', 'data'),
    Output('spatial_dataTable_sparkx', 'columns'),
    # Output('clientside_store_table_sparkx', 'data'),
  ],
  [
    Input('spatial_dropdown_stage_sparkx', 'value'),
    State('spatial_dropdown_featureType_sparkx', 'value'),
  ]
)
def choose_sparkx_plotOptions(stage, featureType):
  if not stage:
    raise PreventUpdate
  df = genes_all_pval.loc[stage].reset_index().rename(columns={"index": "gene"})
  df['id'] = df['gene']
  return (df.to_dict('records'), 
          [{"name": i, "id": i, "deletable": False, 'type': 'numeric', 
              'format':Format(precision=2, scheme=Scheme.exponent)} 
             if i != 'gene' else
             {"name": i, "id": i, "deletable": False} 
             for i in df.columns if i != 'id'],
          # {'data': df.to_dict('records'),
          # 'columns': [{"name": i, "id": i, "deletable": False, 'type': 'numeric', 
          #               'format':Format(precision=2, scheme=Scheme.exponent)} 
          #            if i != 'gene' else
          #              {"name": i, "id": i, "deletable": False} 
          #            for i in df.columns if i != 'id']},
          )

@app.callback(
  Output('spatial_dataTable_sparkx', 'page_size'),
  Input('spatial_dropdown_pageSize_sparkx', 'value'),
)
def update_sprkx_tablePageSize(pageSize):
  if not pageSize:
    raise PreventUpdate
  else:
    return pageSize
  
@app.callback(
  [
    Output('spatial_dataTable_sparkx', 'style_data_conditional'),
  ],
  [
    Input('spatial_dataTable_sparkx', 'derived_viewport_selected_row_ids'),
  ],
  prevent_initial_call=True
)
def update_highlightRow_store(rows):
  if not rows:
    raise PreventUpdate
  row_highlighting = [
      {
          'if': {"filter_query": "{{id}} ={}".format(i)},
          'background_color': 'tomato',
          'color': 'white'
      } for i in rows
  ]
  return (row_highlighting,)


@app.callback(
  [
    Output('spatial_featureGraph_sparkx', 'figure'),
  ],
  [
    Input('spatial_dataTable_sparkx', 'selected_row_ids'),
    Input('spatial_dropdown_stage_sparkx', 'value'),
    Input('spatial_dropdown_featureType_sparkx', 'value')
  ],
  prevent_initial_call=True,
  background=True,
  manager=background_callback_manager,
)
def choose_sparkx_featureToPlot(row_id, stage, featureType):
  if not row_id:
    raise PreventUpdate
  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_mtx[stage]

  return (
    show_feature_spatial(adata, row_id[0], embedding = adata.obs[['x_flatten', 'y_flatten']], sort=True, ascending=True,
                            facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']}),
  )

# @app.callback(
#   Output('spatial_dataTable_sparkx', 'style_data_conditional'),
#   Input('spatial_dataTable_sparkx', 'derived_viewport_selected_row_ids'),
# )
# def highlight_dataTable_row(rows):
#   if not rows:
#     raise PreventUpdate
#   row_highlighting = [
#       {
#           'if': {"filter_query": "{{id}} ={}".format(i)},
#           'background_color': 'tomato',
#           'color': 'white'
#       } for i in rows
#   ]
#   return row_highlighting



@app.callback(
  [
    Output('spatial_celltypeGraph_sparkx', 'figure'),
  ],
  [
    Input('spatial_dropdown_stage_sparkx', 'value'),
  ],
  background=True,
  manager=background_callback_manager,
)
def update_sparkx_stageCtpGraph(stage):
  if not stage:
    raise PreventUpdate
  return (
      show_celltype_spatial(exp_data[stage], embedding = exp_data[stage].obs[['x_flatten', 'y_flatten']],
                           facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']}),
  )


# ## all layout

# In[297]:


# controls = dbc.Col(
#   spatial_controls,
#   id = 'controls'
# )
tabs = dbc.Col(
  spatial_tabs,
  id = 'tabs'
)


# In[298]:


app.layout = dbc.Container(
  [
    header,
    dbc.Row([
      # dbc.Col([
      #   controls,
      # ], width=2),
      dbc.Col([
        tabs,
      ], width=12)
    ],)
  ],
  fluid=True,
  className="Container-all",
)


# # run app

# In[299]:


app.config.suppress_callback_exceptions = True
if __name__ == "__main__":
    app.run(
    host='10.193.0.208',
    port='8053',
    debug=True,
    use_reloader=False,
    jupyter_mode = 'external'
)


# In[ ]:




