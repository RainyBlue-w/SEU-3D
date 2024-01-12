#!/usr/bin/env python
# coding: utf-8

import dash
dash.register_page(__name__, path='/')

# In[] env

from math import isnan
import math
from functools import reduce
from dash import Dash, dcc, html, dash_table, no_update, State, Patch, DiskcacheManager, clientside_callback, ctx, ClientsideFunction
from dash import ALL, MATCH, ALLSMALLER
from dash.dash_table.Format import Format, Group, Scheme, Symbol
from dash.exceptions import PreventUpdate
import dash_daq as daq
import dash_ag_grid as dag
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import feffery_antd_components as fac
import feffery_utils_components as fuc
import dash_bootstrap_components as dbc
from dash_extensions.enrich import Output, Input, html, callback, DashProxy, LogTransform, DashLogger, Serverside, ServersideOutputTransform
from dash_extensions.enrich import MultiplexerTransform

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


# In[] data

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

# In[] functions

def show_celltype_spatial(adata, embedding, cmap = ctp_cmap, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  pdf = pd.DataFrame(np.array(embedding), 
                    index=adata.obs_names, 
                    columns=['x', 'y'])
  pdf = pd.concat([pdf, adata.obs[['celltype','germ_layer']]], axis=1)
  pdf = pdf.sort_values(by='celltype')
  plot = px.scatter(
  	data_frame = pdf,
    x = 'x', y = 'y', color = 'celltype',
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

def show_feature_spatial(adata, feature, embedding, cmap = None, sort=False, ascending=True, **kws):
    embedding = embedding.loc[adata.obs_names,:]
    if cmap is None:
        cmap = [(0.00, "rgb(244,244,244)"),
                (0.05, "rgb(244, 244, 244)"),
                (1.00, "rgb(34, 94, 168)")
                ]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['x', 'y'])
    pdf = pd.concat([pdf, adata[:,feature].to_df(), adata.obs['germ_layer']], axis=1)
    if sort is True:
      pdf = pdf.sort_values(by=feature, ascending=ascending)
    plot = px.scatter(
    		data_frame = pdf,
    		x = 'x', y = 'y', color = feature,
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

def show_celltype_spatial_3D(adata, embedding, cmap = ctp_cmap, **kws):
    embedding = embedding.loc[adata.obs_names,:]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['x', 'y', 'z'])
    pdf = pd.concat([pdf, adata.obs['celltype']], axis=1)
    pdf = pdf.sort_values(by='celltype')
    plot = px.scatter_3d(
    		data_frame = pdf,
    		x = 'x', y = 'y', z='z',
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
                      marker_opacity = 0.9)
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

def vector_to_rgba(v):
  color = list(v.keys())
  color = [str(math.ceil(v[i])) if i in color else '255' for i in ['R', 'G', 'B'] ]
  if(all([ i=='255' for i in color])):
    rgba = 'rgba(255,255,255,0)'
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

  delta = 255 - exp.div(exp.max(axis=0), axis=1)*255
  delta[others] = 255

  def color_geoMean(a,b):
    a = 255-a
    b = 255-b
    geoMean = numpy.sqrt((a**2+b**2)/2)
    # geoMean = ((a**3+b**3)/2)**(1/3)
    color = 255 - geoMean
    return color
  def mean(a,b, c=None):
    if c:
      return (a+b+c)/3
    else:
      return (a+b)/2

  if len(colors)==1:
    color = pd.DataFrame({
        colors[0] : 255,
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

def show_multiFeatures_spatial(adata, features_dict, embedding, **kws):
  
  embedding = embedding.loc[adata.obs_names,:]
  pdf = pd.DataFrame(np.array(embedding), 
                      index=adata.obs_names, 
                      columns=['x', 'y'])
  features = list(features_dict.values())
  pdf = pd.concat([pdf, adata[:,features].to_df(), adata.obs['germ_layer']], axis=1)
  colors = multiGenes_show_color(adata, features_dict)
  
  plot = go.Figure()
  plot.add_trace(go.Scatter(
    x = pdf.x, y = pdf.y, mode='markers',
    marker={'color': colors}
  ))
  plot.update_yaxes(visible=False)
  plot.update_xaxes(visible=False)
  plot.update_traces(marker_size=4.5,
                    marker_opacity=1)
  plot.update_layout(
    margin=dict(l=10, r=10, t=30, b=0),
    plot_bgcolor = '#ffffff', 
    uirevision='constant',
  )
  plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                              font_size = 20)) 
  return plot

def show_multiFeatures_spatial_3D(adata, features_dict, embedding, **kws):
  
  embedding = embedding.loc[adata.obs_names,:]
  pdf = pd.DataFrame(np.array(embedding), 
                      index=adata.obs_names, 
                      columns=['x', 'y', 'z'])
  
  features = [i for i in features_dict.values() if i]
  pdf = pd.concat([pdf, adata[:,features].to_df()], axis=1)

  colors = multiGenes_show_color(adata, features_dict)
  plot = go.Figure()
  plot.add_trace(go.Scatter3d(
    x = pdf.x, y = pdf.y, z = pdf.z, 
    mode='markers', marker={'color': colors}
  ))
  plot.update_traces(marker_size=3,
                    marker_opacity = 0.9)
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

def show_features_spatial_regularExp(adata, stage,  odir, featureType, embedding, pattern=None, features=None,cmap=None, sort=False, ascending=True, dpi=100, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  if not features:
    
    img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features = [i for i in features if i in genes_min_pval.loc[stage].index.to_list()]

  else:
    img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)

  ordered_features = genes_min_pval.loc[stage].loc[features].sort_values(by='PadjMinVal', ascending=True).index.to_list()

  pdf = pd.DataFrame(np.array(embedding), 
                      index=adata.obs_names, 
                      columns=['x', 'y'])
  
  features_df = adata[:,ordered_features].to_df()
  features_df.columns = ['%s\n(%.2e)' % (i,genes_min_pval.loc[stage].loc[i, 'PadjMinVal']) for i in features_df.columns]

  pdf = pd.concat([pdf, features_df, adata.obs.germ_layer], axis=1)
  pdf = pd.melt(pdf, id_vars = ['x', 'y', 'germ_layer'], var_name='feature', value_name='value')
  if sort:
    pdf = pdf.sort_values(by='value', ascending=ascending)
  pdf['feature'] = pdf['feature'].astype('category').values.reorder_categories(features_df.columns)
  pdf['germ_layer'] = pdf['germ_layer'].astype('category').values.reorder_categories(['ectoderm','mesoderm','endoderm'])

  # plotnine
  (
    ggplot() + 
    geom_point(data = pdf, mapping=aes(x='x', y='y', color='value')) +
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

def show_featuresCtpcounts_spatial_regularExp(adata, stage, odir, featureType, embedding, pattern=None, features=None, cmap=ctp_cmap, sort=False, ascending=True, dpi=150, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  if not features:
    img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features = [i for i in features if i in genes_min_pval.loc[stage].index.to_list()]

  else:
    img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)

  ordered_features = genes_min_pval.loc[stage].loc[features].sort_values(by='PadjMinVal', ascending=True).index.to_list()
  
  ctp_counts = {}
  for gene in ordered_features:
    df = adata[:,gene].to_df()
    # thr = df.min()+(df.max()-df.min())*0.05
    df = df[df[gene] > 0]
    counts = pd.DataFrame(adata.obs['celltype'][df.index].value_counts())
    counts['gene'] = gene
    counts['count'] = counts['celltype']/sum(counts['celltype'])
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

def show_metadata_discrete_spatial_germLayer(adata, stage,  odir, embedding, obs=['celltype'], cmap = ctp_cmap, dpi=150, sort = False, ascending=True, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  img_dir = '/rad/wuc/dash_data/spatial/tmp/%s/%s/poltFeatureSeries_%s_plot_dpi%d.png' % (odir, stage, obs, dpi)

  # if os.path.exists( img_dir ):
  #   return img_dir

  pdf = pd.DataFrame(np.array(embedding), 
                     index=adata.obs_names, 
                     columns=['x', 'y'])
  
  obs_df = pd.DataFrame(adata.obs[obs])

  pdf = pd.concat([pdf, obs_df, adata.obs.germ_layer], axis=1)
  pdf = pd.melt(pdf, id_vars = ['x', 'y', 'germ_layer'], var_name='obs', value_name='value')
  if sort:
    pdf = pdf.sort_values(by='value', ascending=ascending)
  pdf['obs'] = pdf['obs'].astype('category').values.reorder_categories(obs_df.columns)

  # plotnine
  (
    ggplot() + 
    geom_point(data = pdf, mapping=aes(x='x', y='y', color='value'), size=0.4) +
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

def pre_SubsetCellsByObs(
  obs : pd.DataFrame, 
  germ_layers : List[str],
  celltypes : List[str],
  xyz_range : Dict,
  ):
  '''
  germ_layers: ectoderm, mesoderm, endoderm
  xy_range: dict(x=(minX, maxX),y=(minY, maxY))
  '''
  subset = obs.query(
    '(germ_layer in @germ_layers) & \
    (celltype in @celltypes) & \
    ( @xyz_range["x"][0] <= x <= @xyz_range["x"][1] ) & \
    ( @xyz_range["y"][0] <= y <= @xyz_range["y"][1] ) & \
    ( @xyz_range["z"][0] <= z <= @xyz_range["z"][1] )'
  )
  return subset.index

def show_expViolin(adata, feature, **kws):
  data = adata[:,feature].to_df()[feature]
  # data = data[data>0]
  fig = go.Figure(
    data = go.Violin(
      x=data, y0=f'{feature}({len(data)})', box_visible=True, 
      line_color='black', meanline_visible=True,
      fillcolor='lightseagreen', opacity=0.6,
      jitter=0.1, **kws,  marker_size=2.5, points='outliers',
    )
  )
  
  fig.update_traces(orientation='h', side='positive', width=1.5)
  fig.update_layout(
    plot_bgcolor = 'rgba(200,200,200,0.1)', showlegend=False
  ).update_yaxes(
    gridcolor='rgba(200,200,200,0.6)', gridwidth=1,
  ).update_xaxes(
    dtick=1, gridcolor='#ffffff', gridwidth=1, griddash='solid'
  )
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

def show_multiFeatures_expViolin(adata, features_dict, **kws):
  
  fig = go.Figure()
  
  filt_dict = {}
  for color,feature  in features_dict.items():
      if feature:
          filt_dict[color] = feature

  for color in list(filt_dict.keys())[::-1]:
    data = adata[:,filt_dict[color]].to_df()[filt_dict[color]]
    fig.add_trace(
      go.Violin(
        x=data, y0=f'{filt_dict[color]}({len(data)})', box_visible=False, 
        line_color='black', meanline_visible=False,
        fillcolor=color, opacity=0.6, marker_size=2.5, points='outliers',
         jitter=0.1, **kws
      )
    )
  fig.update_traces(orientation='h', side='positive', width=1.5)
  fig.update_layout(
    plot_bgcolor = 'rgba(200,200,200,0.1)', showlegend=False
  ).update_yaxes(
    gridcolor='rgba(200,200,200,0.6)', gridwidth=1,
  ).update_xaxes(
    dtick=1, gridcolor='#ffffff', gridwidth=1, griddash='solid'
  )
  
  return fig

def show_multiFeatures_ctpExpViolin(adata, features_dict, *kws):
  
  from plotly.subplots import make_subplots
  
  filt_dict = {}
  for color,feature  in features_dict.items():
      if feature:
          filt_dict[color] = feature
  features = list(filt_dict.values())

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

  delta = exp.div(exp.max(axis=0), axis=1)*244
  delta[others] = 0

  def delta_geoMean(a,b):
    geoMean = numpy.sqrt((a**2+b**2)/2)
    # geoMean = ((a**3+b**3)/2)**(1/3)
    return geoMean
  def mean(a,b, c=None):
    if c:
      return (a+b+c)/3
    else:
      return (a+b)/2

  if len(colors)==1:
    color = pd.DataFrame({
        colors[0] : 244,
        others[0] : 244-delta[colors[0]],
        others[1] : 244-delta[colors[0]],
    })
  elif len(colors)==2:
    color = pd.DataFrame({
        colors[0] : 244-delta[colors[1]],
        colors[1] : 244-delta[colors[0]],
        others[0] : 244-delta_geoMean(delta[colors[1]],delta[colors[0]]),
    })
  elif len(colors)==3:
    color = pd.DataFrame({
        'R' : 244-delta_geoMean(delta['G'], delta['B']),
        'G' : 244-delta_geoMean(delta['R'], delta['B']),
        'B' : 244-delta_geoMean(delta['R'], delta['G']),
    })
  
  color['RGBA'] = color.apply(vector_to_rgba, axis=1)
  return color['RGBA']

def hex_to_rgbList(hex_color):
  hex_color = hex_color.replace(' ', '').replace('#', '')
  if len(hex_color) == 6:
    r = int(hex_color[0:2], 16)
    g = int(hex_color[2:4], 16)
    b = int(hex_color[4:6], 16)
  return [r,g,b]

def mix_multipy(color, alpha):

  def multipy(x,y):
    return x*y/255

  def mix(x, y):
    alpha = x[3]+y[3]-x[3]*y[3]
    if alpha==0:
      return [244,244,244, 0]
    else:
      R = np.round( (x[3]*(1-y[3])*x[0]+x[3]*y[3]*multipy(x[0],y[0])+(1-x[3])*y[3]*y[0])/alpha).astype(int)
      G = np.round( (x[3]*(1-y[3])*x[1]+x[3]*y[3]*multipy(x[1],y[1])+(1-x[3])*y[3]*y[1])/alpha).astype(int) 
      B = np.round( (x[3]*(1-y[3])*x[2]+x[3]*y[3]*multipy(x[2],y[2])+(1-x[3])*y[3]*y[2])/alpha).astype(int)
      return [R,G,B,alpha]

  array = []
  for c,a in zip(color, alpha):
    array.append(c.copy())
    array[-1].append(a)

  res = reduce(mix, array)
  res = f'rgb{res[0],res[1],res[2]}'

  return res

def color_mixer(adata, genes_dict):
  import numpy
  genes_dict_copy = genes_dict.copy()
  _ = [genes_dict_copy.pop(color) for color in genes_dict.keys() if not genes_dict[color]]
  colors = [hex_to_rgbList(c) for c in genes_dict_copy.keys()]
  genes = list(genes_dict_copy.values())
  
  exp = adata[:,genes].to_df()
  
  alpha = exp.div(exp.max(axis=0), axis=1)
  
  cell_colors = alpha.apply( axis=1, func=lambda row: mix_multipy(colors,row))
  
  return cell_colors

def cal_moran_3D(adata):
  tmp = adata.copy()
  sq.gr.spatial_neighbors(tmp, spatial_key='X_spatial')
  sq.gr.spatial_autocorr(tmp, mode='moran', n_jobs=1)
  df = tmp.uns['moranI'][['I']]
  df.columns = ["Moran's I"]
  return df

# In[] global vars
colorPicker_swatches = [
  "#25262b", "#868e96", "#fa5252", "#e64980", "#be4bdb", "#7950f2", "#4c6ef5",
  '#225ea8', "#228be6", "#15aabf", "#12b886", "#40c057", "#82c91e", "#fab005", "#fd7e14",
]

initColor_multiName = [
  "#fa5252", "#228be6", "#40c057", "#fd7e14", "#be4bdb", "#e64980", "#15aabf", "#fab005", "#868e96", 
]


# In[] app/wdigets/flat :

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
            clearable=True,
            searchable=True,
        ),
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


# In[] app/widgets/3D:

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
    dcc.Store(id='STORE_cellsCtpFilter_3D'),
    dcc.Store(id='STORE_cellsIntersection_3D'),
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

def iconHover_colorPicker(init_color: str, id: Dict, swatches: List[str], placement='left', trigger='click'):
  return fac.AntdPopover(
    # openDelay=200,
    placement = placement,
    trigger= trigger,
    children=[
      dmc.ActionIcon(
        DashIconify(icon = 'fluent:circle-48-filled', color=init_color, width=48),
        variant='transparent', id=id['dmc_ActionIcon'], mt=3
      ),
    ],
    content = [
      dmc.ColorPicker(id=id['dmc_ColorPicker'], format='hex', value=init_color, swatches=swatches),
      dmc.TextInput(value=init_color, id=id['dmc_TextInput']),
    ]
  )

def drawerContent_ctpColorPicker(celltypes: List[str], cmap: Dict, swatches=colorPicker_swatches):
  stack = dmc.Stack(
    children=[
      dmc.Grid(
        gutter = 2,
        children=[
          dmc.Col(dmc.Text(ctp), span=10),
          dmc.Col(
            iconHover_colorPicker(
              id = {
                'dmc_ActionIcon': {'type': 'ACTIONICON_colorCtp_3D', 'id': ctp},
                'dmc_ColorPicker': {'type': 'COLORPICKER_colorCtp_3D', 'id': ctp},
                'dmc_TextInput': {'type': 'TEXT_colorCtp_3D', 'id': ctp},
              },  placement='right', init_color=cmap[ctp], swatches = [ctp_cmap[ctp]]+swatches
            ),
            span=2
          )
        ],
      )
      for ctp in celltypes
    ],
  )

  return stack

# In[] app/widgets/series:

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
spatial_textarea_featureLists_series = html.Div(
  [
    dbc.Label("name list:"),
    dbc.Col(
      [
        dbc.Textarea(id = "spatial_textarea_featureLists_series",
          placeholder="paste feature names here(seperated by any charactor)\ne.g.  A  B,C,  D\nE##G_H,#I@KK%%G%(U),(V)|(W)\"X\",\"Y\"^Q*I",
          rows=8, className="mb-3",),
        dbc.Button('Plot', id='spatial_inputButton_featureLists_series_plot',
                          n_clicks=0, color='primary'),
      ]
    )
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
# In[] app/tabs/:

spatial_tab_plotFeature = dbc.Tab(
  [dbc.Row([
    dbc.Col([
      dbc.Card(
        [
          spatial_dropdown_featureType,
          spatial_dropdown_featureName,
          # spatial_dropdown_plotDim,
          spatial_dropdown_stage
        ],
        body=True,
        id = 'spatial_control_plotFeature'
      ),      
    ], width=2),
    dbc.Col([dbc.Row([
      dbc.Col([
        dcc.Graph(figure={}, id="spatial_plotFeature_graph", style={'height': "40vh"}),
      ],align = "center",className="g-0", width=12),
      dbc.Col([
        dcc.Graph(figure={}, id="spatial_plotFeature_graph_ctp", style={'height': "40vh"}),
      ],align = "center",className="g-0", width=12)
    ]),],width=10),
  ])],
  label = "Plot feature(flat)",
  tab_id = "spatial_tab_plotFeature"
)

spatial_tab_plotFeature3D = dbc.Tab(
  [dmc.Grid([
    # options
    dmc.Col([
      fuc.FefferySticky([
        fac.AntdSpace(
          size=0,
          direction='vertical',
          className='fac-AntdSpace-sideBar',
          children=[
            # Select data
            fac.AntdCollapse(
              isOpen = True,
              forceRender = True,
              className = 'fac-AntdCollapse-sidebar',
              ghost=True,
              title = dmc.Text('Select data', className='dmc-Text-sidebar-title', style={'font-size': 18}),
              children = [
                dmc.Grid([
                  dmc.Col([
                    dbc.Label("Feature type"),
                    dcc.Dropdown(
                      ['Gene', 'Regulon'],
                      'Gene',
                      id="DROPDOWN_featureType_3D",
                      clearable=False,
                      searchable=True,
                    ),
                  ], span=6),
                  dmc.Col([
                    dbc.Label("Stage"),
                    dcc.Dropdown(
                      ['E7.5', 'E7.75', 'E8.0'],
                      'E7.75',
                      id="DROPDOWN_stage_3D",
                      clearable=False,
                      searchable=True,
                    ),
                  ], span=6),
                  dmc.Col(dmc.Text(id='TEXT_dataSummary_3D', color='gray'), span=12),
                  fac.AntdCollapse(
                    title = 'germ layer', className='fac-AntdCollapse-inline',
                    forceRender=True, isOpen=True, ghost=True,
                    children = [
                      dmc.ChipGroup(
                        children = [
                          dmc.Chip(
                            x.capitalize(),
                            value=x
                          ) for x in ['ectoderm', 'mesoderm', 'endoderm']
                        ],
                        value = ['ectoderm', 'mesoderm', 'endoderm'],
                        id = 'CHIPGROUP_germLayer_3D',
                        multiple = True,
                      )
                    ]
                  ),
                ], gutter='xs'),
              ]
            ),
            # Plot options
            fac.AntdCollapse(
              isOpen = True,
              forceRender = True,
              className = 'fac-AntdCollapse-sidebar',
              ghost=True,
              title = dmc.Text('Plot options', className='dmc-Text-sidebar-title', style={'font-size': 18}),
              children = [          
                dmc.Tabs(
                  [
                    dmc.TabsList([
                      dmc.Tab('Settings', value='settings'),
                      dmc.Tab('Single', value='single'),
                      dmc.Tab('Multiple', value='multi'),
                    ], grow=False),
                    # settings
                    dmc.TabsPanel(
                      [
                        dmc.Text('Point size', className='dmc-Text-label'),
                        dmc.NumberInput(value=3, step=0.5, min=0.1, id='NUMBERINPUT_pointSize_3D', precision=1),
                        dmc.Space(h=10),
                        dmc.Switch(label='Hide axes', id='SWITCH_hideAxes_3D', size='md',
                                  onLabel=DashIconify(icon='solar:eye-closed-linear', width=14), 
                                  offLabel=DashIconify(icon='solar:eye-linear', width=14) ),
                        dmc.Space(h=10),
                        dmc.Switch(label='Hide non-expressing cells', id='SWITCH_hideZero_3D',  size='md',
                                  onLabel=DashIconify(icon='solar:eye-closed-linear', width=14), 
                                  offLabel=DashIconify(icon='solar:eye-linear', width=14) ),
                        dmc.Space(h=10),
                        dmc.Text('Projection type', className='dmc-Text-label'),
                        dmc.SegmentedControl(
                          value='perspective', 
                          data=[
                            {'value': 'perspective', 'label': 'Perspective'},
                            {'value': 'orthographic', 'label': 'Orthographic'},
                          ], 
                          fullWidth=True, id='SEGMENTEDCONTROL_projection_3D'
                        ),
                      ],
                      value = 'settings',
                    ),
                    # single
                    dmc.TabsPanel(
                      [
                        dmc.Grid([
                          dmc.Col([
                            dcc.Dropdown(
                              options = exp_data['E7.75'].var_names,
                              value = 'Cdx1',
                              id="DROPDOWN_singleName_3D",
                              clearable=False
                            ),
                          ], span=10),
                          dmc.Col([
                            iconHover_colorPicker(
                              id={
                                'dmc_ActionIcon': 'ACTIONICON_colorSingle_3D', 
                                'dmc_ColorPicker': 'COLORPICKER_single_3D', 
                                'dmc_TextInput': 'TEXT_colorSingle_3D',
                                },
                              init_color='#225ea8', swatches=colorPicker_swatches,
                            )
                          ], span=2),
                          dmc.Col([
                            dmc.Button('Plot', id='BUTTON_singlePlot_3D', color='dark', fullWidth=True,
                                      leftIcon=DashIconify(icon="gis:cube-3d", width=24)),
                          ], span=12),
                        ], gutter='xs')  
                      ],
                      value='single',
                    ),
                    # multi
                    dmc.TabsPanel(
                      [
                        # extendable selector
                        html.Div(
                          [
                            dmc.Grid([
                              dmc.Col(dcc.Dropdown(options = [], id={'type': 'DROPDOWN_multiName_3D', 'index': 0}), span=10),
                              dmc.Col(
                                iconHover_colorPicker(
                                  id={
                                    'dmc_ActionIcon': {'type':'ACTIONICON_colorMulti_3D', 'index': 0}, 
                                    'dmc_ColorPicker': {'type': 'COLORPICKER_multi_3D', 'index': 0}, 
                                    'dmc_TextInput': {'type': 'TEXT_colorMulti_3D', 'index': 0},
                                    },
                                  init_color=initColor_multiName[0], swatches=colorPicker_swatches,
                                ),
                                span=2
                              ),
                            ]),
                            dmc.Grid([
                              dmc.Col(dcc.Dropdown(options = [], id={'type': 'DROPDOWN_multiName_3D', 'index': 1}), span=10),
                              dmc.Col(
                                iconHover_colorPicker(
                                  id={
                                    'dmc_ActionIcon': {'type':'ACTIONICON_colorMulti_3D', 'index': 1}, 
                                    'dmc_ColorPicker': {'type': 'COLORPICKER_multi_3D', 'index': 1}, 
                                    'dmc_TextInput': {'type': 'TEXT_colorMulti_3D', 'index': 1},
                                    },
                                  init_color=initColor_multiName[1], swatches=colorPicker_swatches,
                                ),
                                span=2
                              ),
                            ]),
                          ],
                          id = 'DIV_multiNameDynamic_3D',
                        ),
                        dcc.Store(data=2, id='STORE_multiNameCurNumber'),
                        # buttons
                        dmc.Grid(
                          [
                            dmc.Col(dmc.Button(
                              'Add', id='BUTTON_addFeature_3D', color='teal', fullWidth=True,
                              leftIcon=DashIconify(icon="fluent:add-square-20-regular", width=20)
                            ), span=23),
                            dmc.Col(dmc.Button(
                              'Delete', id='BUTTON_deleteFeature_3D', color='red', fullWidth=True,
                              leftIcon=DashIconify(icon="fluent:subtract-square-20-regular", width=20)
                            ), span=27),
                            dmc.Col(dmc.Button(
                              'Plot', id='BUTTON_multiPlot_3D', color='dark', fullWidth=True,
                              leftIcon=DashIconify(icon="gis:cube-3d", width=24),
                            ), span=50),
                          ],
                          columns=50,
                        ),
                        dcc.Store(id='STORE_multiNameInfo_3D'),
                      ],
                      value='multi',
                    ),
                  ], 
                  orientation = 'vertical',
                  placement = 'left',
                  className = 'dmc-Tabs-vertical',
                  value = 'single'
                ),
              ]
            ),
            # Slicer
            fac.AntdCollapse(
              isOpen = False,
              forceRender = True,
              className = 'fac-AntdCollapse-sidebar',
              ghost=True,
              title = dmc.Text('Slicer', className='dmc-Text-sidebar-title', style={'font-size': 18}),
              children = [
                dmc.Grid([
                  dmc.Col(dmc.Text('x', className='.dmc-Text-label-center'), span=2),
                  dmc.Col(
                    dcc.RangeSlider(
                      -400,480,10, value=(-400,480), id='SLIDER_Xrange_3D',
                      marks=None, tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    span=10
                  ),
                  dmc.Col(dmc.Text('y', className='.dmc-Text-label-center'), span=2),
                  dmc.Col(
                    dcc.RangeSlider(
                      -410, 470,10, value=(-410,470), id='SLIDER_Yrange_3D',
                      marks=None, tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    span=10
                  ),
                  dmc.Col(dmc.Text('z', className='.dmc-Text-label-center'), span=2),
                  dmc.Col(
                    dcc.RangeSlider(
                      0, 700,10, value=(0,700), id='SLIDER_Zrange_3D',
                      marks=None, tooltip={'placement': 'bottom', 'always_visible': True}
                    ),
                    span=10
                  ),
                  dmc.Col(
                    dmc.Switch(size='md', id='SWITCH_previewBox_3D', label='Preview', checked=False),
                    span=12
                  ),
                  dmc.Col(
                    dmc.Button(
                      'Slice', color='red', id='BUTTON_slice_3D', fullWidth=True, 
                      leftIcon=DashIconify(icon='fluent:screen-cut-20-regular', width=20),
                    ),
                    span=6
                  ),
                  dmc.Col(
                    dmc.Button(
                      'Recover', color='teal', id='BUTTON_recover_3D', fullWidth=True, 
                      leftIcon=DashIconify(icon='fluent:arrow-sync-circle-20-regular', width=20),
                    ),
                    span=6
                  )
                ]),
                SET_STORE_Ranges_3D,
              ],
            ),
            # Moran
            fac.AntdCollapse(
              isOpen = False,
              forceRender = True,
              className = 'fac-AntdCollapse-sidebar',
              ghost=True,
              title = dmc.Text('Compute SVG(moran)', className='dmc-Text-sidebar-title', style={'font-size': 18}),
              children = [
                dmc.Grid([
                  dmc.Col(
                    dmc.Button('Compute', id='BUTTON_calMoran_3D', color='dark', fullWidth=True,
                        leftIcon = DashIconify(icon='fluent:clipboard-math-formula-20-regular', width=20) ),
                    span=6,
                  ),
                  dmc.Col(
                    dmc.Button('Result', id='BUTTON_showMoran_3D', fullWidth=True,
                        leftIcon = DashIconify(icon='fluent:clipboard-checkmark-20-regular', width=20) ),
                    span=6,
                  ),
                  dmc.Text('Using current cells to compute SVGs', color='gray', italic=True),
                ]),
                dbc.Offcanvas(
                  [dash_table.DataTable(
                    id='DATATABLE_moranRes_3D',
                    sort_action="native", page_action='native', filter_action="native",
                    page_current= 0, page_size= 20, fill_width=True,
                    style_cell={'textAlign': 'center'},
                    style_table={'overflowX': 'auto'},
                  )],
                  title = 'SVGs:',
                  placement='end', scrollable=True, backdrop=False, is_open=False,
                  id = 'OFFCANVAS_moranRes_3D',
                ),
              ],
            ),
          ], 
        ),
      ], top=10),
    ], span=9),
    # viewer
    dmc.Col([
      SET_STORE_JSONtoPlot_3D,
      # scatter3d
      dmc.Grid([
        dmc.Col([
          dcc.Graph(figure={}, id="FIGURE_3Dexpression", 
                    className='dcc-figure'),
        ], span=20),
        dmc.Col([
          dcc.Graph(figure={}, id="FIGURE_3Dcelltype",
                    className='dcc-figure'),
        ], span=20),
        dmc.Col([
          # DIY-legend
          dmc.Grid([
            # set colors
            dmc.Col(dmc.Button(
              'Setting Colors', variant="gradient", gradient={"from": "grape", "to": "pink", "deg": 35},
              id='BUTTON_setColors_3D', fullWidth=True,
              leftIcon=DashIconify(icon='fluent:color-20-regular', width=20)
            ), span=12),
            # invert selection
            dmc.Col(dmc.Button(
              DashIconify(icon='system-uicons:reverse', width=21), variant='light', color='gray',
              id='BUTTON_invertSelectionCtp_3D', fullWidth=True,
            ), span=4),
            # clear selection
            dmc.Col(dmc.Button(
              DashIconify(icon='fluent:border-none-20-regular', width=20), variant="light", color='gray',
              id='BUTTON_clearSelectionCtp_3D', fullWidth=True,
            ), span=4),
            # all selection
            dmc.Col(dmc.Button(
              DashIconify(icon='fluent:checkbox-indeterminate-20-regular', width=20), variant="light", color='gray',
              id='BUTTON_allSelectionCtp_3D', fullWidth=True,
            ), span=4),
          ], gutter=2),
          dmc.ChipGroup(
            children=[], value=[], multiple=True, align='center', spacing=1, 
            id = 'CHIPGROUP_celltype_3D', className='dmc-ChipGroup-legend'
          ),
          dcc.Store(id='STORE_allCelltypes_3D'),
          fac.AntdDrawer(children=[],title='Setting colors',id='DRAWER_setColorCtp_3D')
        ], span=10)
      ], columns=50),
      # violin
      dbc.Row([
        dbc.Col([
          dbc.Label( 'Normalized expression in all celltypes(left)'),
          dbc.Label('and in each celltype(right):'),
          dcc.Graph(figure={}, id="FIGURE_expViolin_3D")
        ], align='center', width=4),
        dbc.Col([
          dcc.Graph(figure={}, id="FIGURE_ctpViolin_3D")
        ], align='center', width=8)
      ], style={'overflow-y': 'auto'}, id='test_sticky')
    ],span=41),
  ], columns=50)],
  label = "Plot feature(3D)",
  tab_id = "spatial_tab_plotFeature3D",
)

spatial_tab_plotFeatureSeries = dbc.Tab(
  [dbc.Row([
    dbc.Col([
      dbc.Card(
        [
          dbc.Col([
            spatial_dropdown_featureType_series,
            spatial_dropdown_stage_series,
            spatial_input_featureName_series,
            spatial_textarea_featureLists_series,
          ])
        ],
        body=True,
        style = {'position': 'fixed'},
        id = 'spatial_control_plotFeatureSeries'
      )
    ], width=2),
    dbc.Col([
      dbc.Row(
        [
          dbc.Col(width=2),
          dbc.Col([
            dcc.Graph(figure={}, id="spatial_plotFeatureSeries_graph_ctp", style={'height': "50vh"}),
          ],align = "center",className="g-0", width=8),
        ]
      ),
      dbc.Row([
        dbc.Col(
          [
            html.Img(
              # src = Image.open('/rad/wuc/dash_data/spatial/tmp/plotFeatureSeries/E7.75/plotFeatureSeries_E7.75_^Cdx_Gene_dpi100.png'),
                     id = 'spatial_plotFeatureSeries_img', style = {'width': '90vh'})
          ],
          align = "center", className="g-0", width=7),
        dbc.Col(
          [
            html.Img(
              # src = Image.open('/rad/wuc/dash_data/spatial/tmp/plotFeatureSeries/E7.75/plotFeatureSeries_E7.75_ctpCounts_^Cdx_Gene_dpi150.png'),
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
        html.Div(
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
        dbc.Col([
          dcc.Graph(figure={}, id="spatial_featureGraph_sparkx", style={'height': "36vh", 'width': '130vh'}),
          dcc.Graph(figure={}, id="spatial_celltypeGraph_sparkx", style={'height': "36vh", 'width': '130vh'}),
        ], style = {'position': 'fixed'}),
      ],width=8)
    ]),
  ],
  label = "SPARK-X",
  tab_id = "spatial_tab_plotFeatureSparkx"
)

spatial_tabs = dbc.Card(
    dbc.Tabs(
        [spatial_tab_plotFeature, spatial_tab_plotFeature3D, spatial_tab_plotFeatureSeries, spatial_tab_plotFeatureSparkx],
        active_tab = "spatial_tab_plotFeature",  
        id = "spatial_tabs",
    ),
)

# In[] callbacks/flat :

@callback(
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

@callback(
  Output('spatial_plotFeature_graph', 'figure'),
  State('spatial_dropdown_featureType', 'value'),
  Input('spatial_dropdown_featureName', 'value'),
  Input('spatial_dropdown_stage', 'value'),
  # background=True,
  # manager=background_callback_manager,
)
def update_spatial_plotFeature_graph(featureType, name, stage):
  if name is None:
    raise PreventUpdate
  
  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]
  else:
    raise PreventUpdate

  return show_feature_spatial(adata, name, embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True,
                            facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']})

@callback(
  Output('spatial_plotFeature_graph_ctp', 'figure'),
  State('spatial_dropdown_featureType', 'value'),
  Input('spatial_dropdown_stage', 'value'),
  background=True,
  manager=background_callback_manager,
)
def update_spatial_plotFeature_ctpGraph(featureType, stage):

  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]
  else:
    raise PreventUpdate

  return show_celltype_spatial(adata, embedding = coord_data[stage][['x_flatten', 'y_flatten']],
                               facet_col = 'germ_layer', 
                               category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']})


# In[] callbacks/3D(new):

# update_dataSummary
@callback(
  Output('TEXT_dataSummary_3D', 'children'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value')
)
def update_dataSummary_3D(featureType, stage):
  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_data[stage]
    
  str = f'{adata.shape[0]}(cells) × {adata.shape[1]}(features)'
  return str

# update_nameOptions
@callback(
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

@callback(
  Output({'type': 'DROPDOWN_multiName_3D', 'index': MATCH}, 'options'),
  Input({'type': 'DROPDOWN_multiName_3D', 'index': MATCH}, 'search_value'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value'),
  prevent_initial_call=True,
)
def update_nameOptions_multi_3D(search, featureType, stage):
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

# add & delte components for multiName
@callback(
  Output('DIV_multiNameDynamic_3D', 'children'),
  Output('STORE_multiNameCurNumber', 'data'),
  Input('BUTTON_addFeature_3D', 'n_clicks'),
  Input('BUTTON_deleteFeature_3D', 'n_clicks'),
  State('STORE_multiNameCurNumber', 'data'),
  State('DROPDOWN_featureType_3D', 'value'),
  State('DROPDOWN_stage_3D', 'value'),
  prevent_initial_call = True,
)
def add_components_multiName_3D(add, delete, curNumber, featureType, stage):
  
  if featureType == 'Gene':
    opts = exp_data[stage].var_names
  elif featureType == 'Regulon':
    opts = auc_data[stage].var_names
  
  id = ctx.triggered_id

  nextIndex = curNumber
  nextColor = initColor_multiName[int(nextIndex % len(initColor_multiName))]
  
  patch_children = Patch()
  if 'BUTTON_addFeature_3D' in id:
    patch_children.append(
      dmc.Grid([
        dmc.Col(dcc.Dropdown(options = [], id={'type': 'DROPDOWN_multiName_3D', 'index': nextIndex}), span=10),
        dmc.Col(
          iconHover_colorPicker(
            id={
              'dmc_ActionIcon': {'type':'ACTIONICON_colorMulti_3D', 'index': nextIndex}, 
              'dmc_ColorPicker': {'type': 'COLORPICKER_multi_3D', 'index': nextIndex}, 
              'dmc_TextInput': {'type': 'TEXT_colorMulti_3D', 'index': nextIndex},
              },
            init_color=nextColor, swatches=colorPicker_swatches,
          ),
          span=2
        ),
      ])
    )
    nextNumber = curNumber+1
  elif 'BUTTON_deleteFeature_3D' in id:
    del patch_children[nextIndex-1]
    nextNumber = curNumber-1 if curNumber>0 else 0

  return patch_children, nextNumber

# store_previewRange
clientside_callback(
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
clientside_callback(
  ClientsideFunction(
    namespace='plotFunc_3Dtab',
    function_name='store_sliceRange'),
  Output('STORE_sliceRange_3D', 'data'),
  Input('BUTTON_slice_3D', 'n_clicks'),
  Input('BUTTON_recover_3D', 'n_clicks'),
  State('STORE_previewRange_3D', 'data'),
  State('STORE_maxRange_3D', 'data'),
)

# store_cellsInfo_forJSONtoPlot (download: ~2.5M,500ms ; compute 250ms)

@callback(
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

@callback(
  Output('STORE_cellsObsFilter_3D', 'data'),
  
  Input('STORE_sliceRange_3D', 'data'),
  Input('CHIPGROUP_germLayer_3D', 'value'),
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
  return Serverside(obsnames_filt)

# store_expInfo_forJSONtoPlot (download: <0.43M,<80ms; compute 320ms)
@callback(
  Output('STORE_cellsExpFilter_3D', 'data'),
  Output('STORE_singleExp_3D', 'data'),
  Output('STORE_multiExp_3D', 'data'),
  Output('STORE_ifmulti_3D', 'data'),
  Output('STORE_mixedColor_3D', 'data'),
  
  Input('BUTTON_singlePlot_3D', 'n_clicks'),
  Input('BUTTON_multiPlot_3D', 'n_clicks'),
  Input('SWITCH_hideZero_3D', 'checked'),
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
      cellsExpFilter = exp[(exp>0)[sname]].index.to_list()
    else:
      cellsExpFilter = exp.index.to_list()
    return (ifmulti, exp, cellsExpFilter)
  
  def return_multi():
    ifmulti = True
    mixColor = color_mixer(adata, minfo)
    if hideZero:
      cellsExpFilter = mixColor[mixColor!='rgb(244, 244, 244)'].index.to_list()
    else:
      cellsExpFilter = mixColor.index.to_list()
    return (ifmulti, [], cellsExpFilter, mixColor.to_dict()) 
  
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
        return (Serverside(cellsExpFilter), exp, no_update, ifmulti, no_update)
      else:
        ifmulti,_,cellsExpFilter,mixcolor = return_multi()
        return (Serverside(cellsExpFilter), no_update, no_update, ifmulti, mixcolor)

    elif 'BUTTON_singlePlot_3D' in btn_id:
      ifmulti,exp,cellsExpFilter = return_single()
      exp = exp.to_dict('index')
      if hideZero:
        return (Serverside(cellsExpFilter), exp, no_update, ifmulti, no_update)
      else:
        return (no_update, exp, no_update, ifmulti, no_update)
    
    elif 'BUTTON_multiPlot_3D' in btn_id:
      ifmulti,_,cellsExpFilter,mixcolor = return_multi()
      # exp = return_multiExp()
      if hideZero:
        return (Serverside(cellsExpFilter), no_update, no_update, ifmulti, mixcolor)
      else:
        return (no_update, no_update, no_update, ifmulti, mixcolor)
    
    elif 'SWITCH_hideZero_3D' in btn_id:
      
      if not hideZero:
        cellsExpFilter = adata.obs_names.to_list()
        return (Serverside(cellsExpFilter), no_update, no_update, no_update, no_update)
      
      else:
        if not ifmulti:
          _,_,cellsExpFilter = return_single()
          return (Serverside(cellsExpFilter), no_update, no_update, no_update, no_update)
        else:
          _,_,cellsExpFilter,_ = return_multi()
          return (Serverside(cellsExpFilter), no_update, no_update, no_update, no_update)

  else:
      ifmulti,exp,cellsExpFilter = return_single()
      exp = exp.to_dict('index')
      return (Serverside(cellsExpFilter), exp, no_update, ifmulti, no_update)

# update ChipGroup-celltype chips
@callback(
  Output('CHIPGROUP_celltype_3D', 'children'),
  Output('CHIPGROUP_celltype_3D', 'value'),
  Output('STORE_allCelltypes_3D', 'data'),
  Input('STORE_cellsObsFilter_3D', 'data'),
  Input('STORE_cellsExpFilter_3D', 'data'),
  State('DROPDOWN_stage_3D', 'value'),
  State('STORE_ctpCmap_3D', 'data'),
)
def update_chipGroupCelltype_3D(obsFilter, expFilter, stage, cmap):
  cells = list(set(obsFilter)&set(expFilter))
  cells.sort()
  celltypes = list(exp_data[stage].obs['celltype'][cells].unique())
  
  chips = [
    dmc.Chip(
      children=ctp, value=ctp, size='xs', color='gray', variant='filled', type='radio',
      styles = {
        'label': {
          'color': cmap[ctp],
          'font-size': 12,
          'font-weight': 'bold'
        },
        'checkIcon': {
          'color': cmap[ctp],
        }
      },
      id = {'type': 'CHIP_ctpColorLegend_3D', 'id': ctp}
    ) 
    for ctp in celltypes
  ]
  
  return chips, celltypes, celltypes

@callback(
  Output('CHIPGROUP_celltype_3D', 'value'),
  Input('BUTTON_invertSelectionCtp_3D', 'n_clicks'),
  State('CHIPGROUP_celltype_3D', 'value'),
  State('STORE_allCelltypes_3D', 'data'),
  prevent_initial_call=True,
)
def invertSelection_celltypesButton_3D(click, curValue, allCelltypes):
  return list(set(allCelltypes) - set(curValue))

@callback(
  Output('CHIPGROUP_celltype_3D', 'value'),
  Input('BUTTON_clearSelectionCtp_3D', 'n_clicks'),
  prevent_initial_call=True
)
def clearSelection_celltypesButton_3D(click):
  return []

@callback(
  Output('CHIPGROUP_celltype_3D', 'value'),
  Input('BUTTON_allSelectionCtp_3D', 'n_clicks'),
  State('CHIPGROUP_celltype_3D', 'value'),
  State('STORE_allCelltypes_3D', 'data'),
  prevent_initial_call=True
)
def allSelection_celltypesButton_3D(click, curValue, allCelltypes):
  if set(curValue) == set(allCelltypes):
    return no_update
  else:
    return list(set(allCelltypes))

# store_ctpInfo_forJSONtoPLot
@callback(
  Output('STORE_cellsCtpFilter_3D', 'data'),
  Input('CHIPGROUP_celltype_3D', 'value'),
  State('DROPDOWN_stage_3D', 'value')
)
def store_ctpInfo_forJSONtoPlot_3D(selectedCtps, stage):
    
  series = exp_data[stage].obs['celltype']
  series = series[series.isin(selectedCtps)]
  
  return series.index.to_list()

# plot_3Dfigure_exp
clientside_callback(
  ClientsideFunction(
    namespace='plotFunc_3Dtab',
    function_name='exp_3Dscatter',
  ),
  Output("FIGURE_3Dexpression", "figure"),
  Input('STORE_obs_3D', 'data'),
  Input('STORE_cellsIntersection_3D', 'data'),
  Input('STORE_singleExp_3D', 'data'),
  Input('STORE_ifmulti_3D', 'data'),
  Input('STORE_mixedColor_3D', 'data'),
  State('SWITCH_hideAxes_3D', 'checked'),
  State('SWITCH_previewBox_3D', 'checked'),
  State('STORE_previewRange_3D', 'data'),
  State('SEGMENTEDCONTROL_projection_3D', 'value'),
)

# colorpicker for singleExp
@callback(
  Output('ACTIONICON_colorSingle_3D', 'children'),
  Output('FIGURE_3Dexpression', 'figure'),
  Input('COLORPICKER_single_3D', 'value'),
)
def colorpicker_for_singleExp_3D(color):
  patch = Patch()
  patch['layout']['coloraxis']['colorscale'][1][1] = color
  icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
  return icon, patch

@callback(
  Output('TEXT_colorSingle_3D', 'value'),
  Output('COLORPICKER_single_3D', 'value'),
  Input('TEXT_colorSingle_3D', 'value'),
  Input('COLORPICKER_single_3D', 'value'),
  prevent_initial_call=True,
)
def linkage_colorPickerAndTextSingle_3D(value1, value2):
  id = ctx.triggered_id
  if id == 'TEXT_colorSingle_3D':
    if((len(value1)==4) or (len(value1)==7)):
      return no_update, value1
    else:
      raise PreventUpdate 
  else:
    return value2, no_update

# colopicker for multiExp

@callback(
  Output('STORE_multiNameInfo_3D', 'data'),
  Input({'type': 'COLORPICKER_multi_3D', 'index': ALL}, 'value'),
  Input({'type': 'DROPDOWN_multiName_3D', 'index': ALL}, 'value'),
)
def store_multiNameInfo_3D(colors, genes):
  return dict(zip(colors, genes))

@callback(
  Output({'type':'ACTIONICON_colorMulti_3D', 'index': MATCH}, 'children'),
  Output({'type': 'TEXT_colorMulti_3D', 'index': MATCH}, 'value'),
  Output({'type': 'COLORPICKER_multi_3D', 'index': MATCH}, 'value'),
  Input({'type': 'TEXT_colorMulti_3D', 'index': MATCH}, 'value'),
  Input({'type': 'COLORPICKER_multi_3D', 'index': MATCH}, 'value'),
  prevent_initial_call=True,
)
def linkage_colorPickerAndTextMulti_3D(value1, value2):
  id = ctx.triggered_id
  
  if id['type'] == 'TEXT_colorSingle_3D':
    if((len(value1)==4) or (len(value1)==7)):
      color = value1
      icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
      return icon, no_update, color
    else:
      raise PreventUpdate   
  else:
    color = value2
    icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
    return icon, color, no_update

# colorpicker for ctpLegend
@callback(
  Output('DRAWER_setColorCtp_3D', 'visible'),
  Input('BUTTON_setColors_3D', 'n_clicks'),
  prevent_initial_call=True,
)
def setCelltypeColorsInDrawer_3D(click):
  return True

@callback(
  Output('DRAWER_setColorCtp_3D', 'children'),
  Input('STORE_allCelltypes_3D', 'data'),
  State('STORE_ctpCmap_3D', 'data'),
  # prevent_initial_call=True, 
)
def generate_drawerLegendContent_3D(curAllCtps, cmap):
  return drawerContent_ctpColorPicker(curAllCtps, cmap)

@callback(
  Output({'type': 'ACTIONICON_colorCtp_3D', 'id': MATCH}, 'children'),
  Output({'type': 'TEXT_colorCtp_3D', 'id': MATCH}, 'value'),
  Output({'type': 'COLORPICKER_colorCtp_3D', 'id': MATCH}, 'value'),
  Input({'type': 'TEXT_colorCtp_3D', 'id': MATCH}, 'value'),
  Input({'type': 'COLORPICKER_colorCtp_3D', 'id': MATCH}, 'value'),
  prevent_initial_call=True,
)
def syncAndReturn_colorValue_3D(text, picker):


  tid = ctx.triggered_id
  celltype = tid['id']
  
  if tid['type'] == 'TEXT_colorCtp_3D':
    if((len(text)==4) or (len(text)==7)):
      color = text  
      icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
      return icon, no_update, color
    else:
      raise PreventUpdate
  else:
    color = picker
    icon = DashIconify(icon = 'fluent:circle-48-filled', color=color, width=48)
    return icon, color, no_update
  
@callback(
  Output('STORE_ctpCmap_3D', 'data'),
  Input({'type': 'TEXT_colorCtp_3D', 'id': ALL}, 'value'),
  prevent_initial_call=True
)
def update_storeCtpCmap_3D(colors):
    triggered = ctx.triggered
    triggered_id = ctx.triggered_id
    if(len(triggered) > 1):
      raise PreventUpdate
    
    color = triggered[0]['value']
    ctp = triggered_id['id']
    print('triggered:', color, 'triggered_id:', ctp)
    patch = Patch()
    patch[ctp] = color
    return patch

@callback(
  Output('FIGURE_3Dcelltype', 'figure'),
  Output('FIGURE_ctpViolin_3D', 'figure'),
  Input('STORE_ctpCmap_3D', 'data'),
  State('STORE_allCelltypes_3D', 'data'),
  prevent_initial_call=True,
)
def update_figureCtpCmap_3D(cmap, curCtps):
  
  patch_fig=Patch()
  for i in range(0, len(curCtps)):
    patch_fig['data'][i]['marker']['color'] =  cmap[curCtps[i]]
  return patch_fig, patch_fig

@callback(
  Output({'type': 'CHIP_ctpColorLegend_3D', 'id': MATCH}, 'styles'),
  Input({'type': 'TEXT_colorCtp_3D', 'id': MATCH}, 'value'),
  prevent_initial_call=True
)
def update_chipColor_3D(color):
  patch = Patch()
  patch['label']['color'] = color
  patch['checkIcon']['color'] = color
  return patch

# plot_3Dfigure_ctp
clientside_callback(
  ClientsideFunction(
    namespace='plotFunc_3Dtab',
    function_name='ctp_3Dscatter',
  ),
  Output("FIGURE_3Dcelltype", "figure"),
  Input('STORE_obs_3D', 'data'),
  Input('STORE_cellsIntersection_3D', 'data'),
  State('SWITCH_hideAxes_3D', 'checked'),
  State('SWITCH_previewBox_3D', 'checked'),
  State('STORE_previewRange_3D', 'data'),
  State('STORE_ctpCmap_3D', 'data'),
  State('SEGMENTEDCONTROL_projection_3D', 'value'),
)

# sync layout between exp and ctp figure
@callback(
  Output("FIGURE_3Dexpression", "figure"),
  Output("FIGURE_3Dcelltype", "figure"),
  Input("FIGURE_3Dexpression", "relayoutData"),
  Input("FIGURE_3Dcelltype", "relayoutData"),
  prevent_initial_call=True,
  background=False,
  # manager=background_callback_manager
)
def update_relayout(expLayout, ctpLayout):
  tid = ctx.triggered_id
  patch = Patch()
  
  if tid == 'FIGURE_3Dexpression':
    if 'scene.camera' in expLayout:
      patch['layout']['scene']['camera'] = expLayout['scene.camera']
    if 'scene.aspectratio' in expLayout:
      patch['layout']['scene']['aspectmode'] = 'manual'
      patch['layout']['scene']['aspectratio'] = expLayout['scene.aspectratio']
      
      
    return no_update, patch

  elif tid == 'FIGURE_3Dcelltype':
    if 'scene.camera' in ctpLayout:
      patch['layout']['scene']['camera'] = ctpLayout['scene.camera']
    if 'scene.aspectratio' in ctpLayout:
      patch['layout']['scene']['aspectmode'] = 'manual'
      patch['layout']['scene']['aspectratio'] = ctpLayout['scene.aspectratio']

    return patch, no_update
  
  else:
    raise PreventUpdate


# update point size
@callback(
  Output('FIGURE_3Dexpression', 'figure'),
  Input('NUMBERINPUT_pointSize_3D', 'value'),
  prevent_initial_call = True
)
def update_expPointSize_3D(size):
  
  patch = Patch()
  patch['data'][0]['marker']['size'] = size
  
  return patch

@callback(
  Output('FIGURE_3Dcelltype', 'figure'),
  Input('NUMBERINPUT_pointSize_3D', 'value'),
  State('STORE_cellsIntersection_3D', 'data'),
  State('DROPDOWN_stage_3D', 'value'),
  prevent_initial_call = True,
)
def update_ctpPointSize_3D(size, cells, stage):
  
  adata = exp_data[stage]
  celltypes = adata.obs.loc[cells, 'celltype'].unique()
  patch = Patch()
  for i in range(0,len(celltypes)):
    patch['data'][i]['marker']['size'] = size
  return patch


# switch projection type
@callback(
  Output('FIGURE_3Dcelltype', 'figure'),
  Output('FIGURE_3Dexpression', 'figure'),
  Input('SEGMENTEDCONTROL_projection_3D', 'value'),
)
def switch_projectionType(type):
  patch=Patch()
  patch['layout']['scene']['camera']['projection'] = {'type': type}
  return patch, patch


# find intersection of 3-filter
@callback(
  Output('STORE_cellsIntersection_3D', 'data'),
  Input('STORE_cellsObsFilter_3D', 'data'),
  Input('STORE_cellsExpFilter_3D', 'data'),
  Input('STORE_cellsCtpFilter_3D', 'data'),
)
def intersection_of_filter(obsFilter, expFilter, ctpFilter):
  tmp = list(set(obsFilter) & set(expFilter) & set(ctpFilter))
  tmp.sort()
  return tmp

# hide axes
@callback(
  Output("FIGURE_3Dexpression", "figure", allow_duplicate=True),
  Output("FIGURE_3Dcelltype", "figure", allow_duplicate=True),
  Input('SWITCH_hideAxes_3D', 'checked'),
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
@callback(
  Output('FIGURE_3Dexpression', 'figure', allow_duplicate=True),
  Output('FIGURE_3Dcelltype', 'figure', allow_duplicate=True),
  Input('SWITCH_previewBox_3D', 'checked'),
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
@callback(
  Output('FIGURE_expViolin_3D', 'figure'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value'),
  Input('STORE_cellsIntersection_3D', 'data'),
  Input('STORE_ifmulti_3D', 'data'),
  Input('BUTTON_singlePlot_3D', 'n_clicks'),
  Input('BUTTON_multiPlot_3D', 'n_clicks'),
  State('DROPDOWN_singleName_3D', 'value'),
  State('STORE_multiNameInfo_3D', 'data'),
  background = True,
  manager = background_callback_manager,
)
def update_spatial_plotFeature3D_expViolin(featureType, stage, cells, ifmulti, splot, mplot, sname, minfo):
  
  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]
  
  adata = adata[cells]

  if not ifmulti:
    fig = show_expViolin(adata, sname)
  else:
    fig = show_multiFeatures_expViolin(adata, minfo)

  return fig

@callback(
  Output('FIGURE_ctpViolin_3D', 'figure'),
  Input('DROPDOWN_featureType_3D', 'value'),
  Input('DROPDOWN_stage_3D', 'value'),
  Input('STORE_cellsIntersection_3D', 'data'),
  Input('STORE_ifmulti_3D', 'data'),
  Input('BUTTON_singlePlot_3D', 'n_clicks'),
  Input('BUTTON_multiPlot_3D', 'n_clicks'),
  State('DROPDOWN_singleName_3D', 'value'),
  State('STORE_multiNameInfo_3D', 'data'),
  background = True,
  manager = background_callback_manager,
)
def update_spatial_plotFeature3D_ctpExpViolin(featureType, stage, cells, ifmulti, splot, mplot, sname, minfo):
  if featureType == 'Gene':
      adata = exp_data[stage]
  elif featureType == 'Regulon':
      adata = auc_data[stage]

  adata = adata[cells]

  if not ifmulti:
    fig = show_ctpExpViolin(adata, sname)
  else:
    fig = show_multiFeatures_ctpExpViolin(adata, minfo)

  return fig

# moran SVG offcanvas
@callback(
  Output('OFFCANVAS_moranRes_3D', 'is_open'),
  Input('BUTTON_showMoran_3D', 'n_clicks'),
  prevent_initial_call = True
)
def show_moranRes_offcanvas(click):
  if click:
    return True
  
@callback(
  Output('DATATABLE_moranRes_3D', 'data'),
  Output('DATATABLE_moranRes_3D', 'columns'),
  
  Input('BUTTON_calMoran_3D', 'n_clicks'),
  State('STORE_cellsIntersection_3D', 'data'),
  State('DROPDOWN_stage_3D', 'value'),
  State('DROPDOWN_featureType_3D', 'value'),
  prevent_initial_call=True,
  background = True,
  manager = background_callback_manager,
  running = [
    (Output('BUTTON_showMoran_3D', 'disabled'), True, False),
    (Output('BUTTON_calMoran_3D', 'children'), '< 1min', 'Compute'),
    (Output('BUTTON_calMoran_3D', 'loading'), True, False),
    (Output('OFFCANVAS_moranRes_3D', 'is_open'), False, True),
  ]
)
def cal_moranRes(click, cells, stage, featureType):
  
  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_data[stage]
  
  df = cal_moran_3D(adata[cells])
  df = df.reset_index(names='Feature')
  return (df.to_dict('records'),
          [
            {"name": i, "id": i, "deletable": False, 'type': 'numeric', 
              'format':Format(precision=4)} 
            for i in df.columns
          ]
        )

# In[] callbacks/series :

@callback(
  Output('spatial_plotFeatureSeries_img', 'src', allow_duplicate=True),
  Output('spatial_plotFeatureSeries_ctpCounts_img', 'src', allow_duplicate=True),
  State('spatial_dropdown_featureType_series', 'value'),
  State('spatial_input_featureName_series', 'value'),
  Input('spatial_inputButton_featureName_series_plot', 'n_clicks'),
  State('spatial_dropdown_stage_series', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('spatial_inputButton_featureLists_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('spatial_inputButton_featureLists_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('spatial_inputButton_featureLists_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('spatial_textarea_featureLists_series', 'disabled', allow_duplicate=True),True, False),
    (Output('spatial_inputButton_featureName_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('spatial_inputButton_featureName_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('spatial_inputButton_featureName_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('spatial_input_featureName_series', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_spatial_plotFeature_graphSeries_pattern(featureType, pattern, click, stage):
  if pattern is None:
    raise PreventUpdate

  if click:
    pattern = '^'+pattern
    if featureType == 'Gene':
      adata = exp_data[stage]
    elif featureType == 'Regulon':
      adata = auc_mtx[stage]

    # img_dir1 = show_features_spatial_regularExp(adata, stage, pattern, 'plotFeatureSeries', featureType, embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True,)
    # img_dir2 = show_featuresCtpcounts_spatial_regularExp(adata, stage, pattern, 'plotFeatureSeries', featureType, embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True,)
    with futures.ThreadPoolExecutor(max_workers=8) as executor:
      t1 = executor.submit(show_features_spatial_regularExp, adata, stage, 'plotFeatureSeries', featureType, 
                              pattern=pattern, embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True)
      t2 = executor.submit(show_featuresCtpcounts_spatial_regularExp, adata, stage,'plotFeatureSeries', featureType, 
                              pattern=pattern, embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True)
      img_dir1 = t1.result()
      img_dir2 = t2.result()
    return Image.open(img_dir1), Image.open(img_dir2)
  else:
    raise PreventUpdate


@callback(
  Output('spatial_plotFeatureSeries_img', 'src', allow_duplicate=True),
  Output('spatial_plotFeatureSeries_ctpCounts_img', 'src', allow_duplicate=True),
  Output('notifications-container', 'children'),
  State('spatial_dropdown_featureType_series', 'value'),
  State('spatial_textarea_featureLists_series', 'value'),
  Input('spatial_inputButton_featureLists_series_plot', 'n_clicks'),
  State('spatial_dropdown_stage_series', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('spatial_inputButton_featureLists_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('spatial_inputButton_featureLists_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('spatial_inputButton_featureLists_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('spatial_textarea_featureLists_series', 'disabled', allow_duplicate=True),True, False),
    (Output('spatial_inputButton_featureName_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('spatial_inputButton_featureName_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('spatial_inputButton_featureName_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('spatial_input_featureName_series', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_spatial_plotFeature_graphSeries_list(featureType, names, click, stage):
  if names is None:
    raise PreventUpdate

  if click:

    names = re.split(", |,| |\n|\'|\"|#|_|%|$|@|\(|\)|\||^|&", names)
    names = [i for i in names if i]
    names = list(set(names))

    if featureType == 'Gene':
      adata = exp_data[stage]
    elif featureType == 'Regulon':
      adata = auc_mtx[stage]
    
    
    names_out = [i for i in names if (i not in adata.var_names) or (i not in genes_min_pval.loc[stage].index)]
    if(names_out):
      note = dmc.Notification(
        title="Features don't exits",
        id = 'series_list_featureNoExit',
        action = 'show',
        message = ','.join(names_out),
        color='orange',
        icon=DashIconify(icon="akar-icons:circle-alert"),
      )
    else:
      note = no_update

    names = list(set(names) - set(names_out))
    
    with futures.ThreadPoolExecutor(max_workers=8) as executor:
      t1 = executor.submit(show_features_spatial_regularExp, adata, stage,  'plotFeatureSeries', featureType,
                            features=names, embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True)
      t2 = executor.submit(show_featuresCtpcounts_spatial_regularExp, adata, stage,'plotFeatureSeries', featureType,
                            features=names, embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True)
      img_dir1 = t1.result()
      img_dir2 = t2.result()
    return Image.open(img_dir1), Image.open(img_dir2), note
  else:
    raise PreventUpdate

@callback(
  Output('spatial_plotFeatureSeries_graph_ctp', 'figure'),
  Input('spatial_dropdown_stage_series', 'value'),
)
def update_spatial_plotFeatureSeries_ctpGraph(stage):
  fig = show_celltype_spatial(exp_data[stage], embedding = coord_data[stage][['x_flatten', 'y_flatten']],
                             facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']})
  fig.update_layout(
    legend = {
      'orientation': 'h'
    }
  )
  return fig


# In[] callbacks/sparkx:

@callback(
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

@callback(
  Output('spatial_dataTable_sparkx', 'page_size'),
  Input('spatial_dropdown_pageSize_sparkx', 'value'),
)
def update_sprkx_tablePageSize(pageSize):
  if not pageSize:
    raise PreventUpdate
  else:
    return pageSize
  
@callback(
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

@callback(
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
    show_feature_spatial(adata, row_id[0], embedding = coord_data[stage][['x_flatten', 'y_flatten']], sort=True, ascending=True,
                            facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']}),
  )

@callback(
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
      show_celltype_spatial(exp_data[stage], embedding = coord_data[stage][['x_flatten', 'y_flatten']],
                           facet_col = 'germ_layer', category_orders={'germ_layer':['ectoderm','mesoderm','endoderm']}),
  )

# In[] app/run:

tabs = dbc.Col(
  spatial_tabs,
  id = 'tabs'
)

layout = dbc.Container(
    [
      html.Div(id='notifications-container'),
      dbc.Row([
        dbc.Col([
          tabs,
        ], width=12)
      ],)
    ],
  fluid=True,
  className="Container-all",
)