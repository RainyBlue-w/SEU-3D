#!/usr/bin/env python
# coding: utf-8

import dash
dash.register_page(__name__)

# In[]: env
from dash import dcc, html, dash_table, Input, Output, callback, no_update, State, Patch, DiskcacheManager, clientside_callback
from dash.exceptions import PreventUpdate
import plotly.express as px
from plotnine import *
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import scanpy as sc
import pandas as pd
import numpy as np
import loompy as lp
import math
import os
from dash_iconify import DashIconify
import matplotlib 
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import re
import diskcache
from dash_extensions.enrich import html
background_callback_manager = DiskcacheManager(diskcache.Cache("/rad/wuc/dash_data/reik/cache"))
matplotlib.use('agg')

# In[]: data

colors = pd.read_csv("/rad/wuc/Reik/celltype_color.csv")
colors['celltype'] = [re.sub(' |/', '_', i) for i in colors['celltype']]
celltype_colors = dict(zip(colors['celltype'], colors['colour']))

stage_list = ['E7.5','E7.75','E8.0','E8.5','E8.75']

umap = pd.read_csv('/rad/wuc/Reik/all_combine_umap.csv', index_col=0)
exp_data = sc.read_h5ad('/rad/wuc/Reik/anndata.h5ad')
exp_data = exp_data[(exp_data.obs['sample'] != 'E8.5_CRISPR_T_WT')&(exp_data.obs['sample'] != 'E8.5_CRISPR_T_KO')]
sc.pp.normalize_total(exp_data, target_sum=1e4)
sc.pp.log1p(exp_data)

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

# In[]: functions/

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

def show_celltype_umap_series(adata, embedding, cmap = celltype_colors, sort=False, **kws):
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
        orientation='h', yanchor='top',xanchor='left', y=1,x=1,
        font = {'size': 12},
        entrywidth=1, entrywidthmode='fraction'
    ),
    uirevision='constant',
    legend_itemsizing = 'constant'
  )
  plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],
                                      font_size = 20)) 
  return plot

def show_features_reik_regularExp(adata, stage, odir, featureType, embedding, pattern=None, features=None, sort=False, ascending=True, dpi=100, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  if not features:
    
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()
    # features = [i for i in features if i in genes_min_pval.loc[stage].index.to_list()]

  else:
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)

  # ordered_features = genes_min_pval.loc[stage].loc[features].sort_values(by='PadjMinVal', ascending=True).index.to_list()

  pdf = pd.DataFrame(np.array(embedding), 
                      index=adata.obs_names, 
                      columns=['x', 'y'])
  
  features_df = adata[:, features].to_df()
  # features_df = adata[:,ordered_features].to_df()
  # features_df.columns = ['%s\n(%.2e)' % (i,genes_min_pval.loc[stage].loc[i, 'PadjMinVal']) for i in features_df.columns]

  # pdf = pd.concat([pdf, features_df, adata.obs.germ_layer], axis=1)
  pdf = pd.concat([pdf, features_df], axis=1)
  # pdf = pd.melt(pdf, id_vars = ['x', 'y', 'germ_layer'], var_name='feature', value_name='value')
  pdf = pd.melt(pdf, id_vars = ['x', 'y'], var_name='feature', value_name='value')
  if sort:
    pdf = pdf.sort_values(by='value', ascending=ascending)
  pdf['feature'] = pdf['feature'].astype('category').values.reorder_categories(features)
  # pdf['germ_layer'] = pdf['germ_layer'].astype('category').values.reorder_categories(['ectoderm','mesoderm','endoderm'])

  # plotnine
  (
    ggplot() + 
    geom_point(data = pdf, mapping=aes(x='x', y='y', color='value')) +
    facet_grid('feature ~ .') + 
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
      strip_text_y = element_text(size=16, face = 'bold', angle=-90)
    )
  ).save(img_dir, width=12, height=8*len(features), dpi=dpi, 
          limitsize=False, verbose=False)
  return img_dir

def show_featuresCtpcounts_reik_regularExp(adata, stage, odir, featureType, embedding, pattern=None, features=None, cmap=celltype_colors, sort=False, ascending=True, dpi=150, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  if not features:
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()
    # features = [i for i in features if i in genes_min_pval.loc[stage].index.to_list()]

  else:
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)

  # ordered_features = genes_min_pval.loc[stage].loc[features].sort_values(by='PadjMinVal', ascending=True).index.to_list()
  
  ctp_counts = {}
  # for gene in ordered_features:
  for gene in features:
    df = adata[:,gene].to_df()
    # thr = df.min()+(df.max()-df.min())*0.05
    df = df[df[gene] > 0]
    counts = pd.DataFrame(adata.obs['celltype'].loc[df.index].value_counts())
    counts['gene'] = gene
    counts['count'] = counts['count']/sum(counts['count'])
    ctp_counts[gene] = counts
  ctp_counts = pd.concat(ctp_counts, axis=0)
  ctp_counts['celltype'] = np.array(ctp_counts.index.to_list())[:,1]
  ctp_counts['text_y'] = ctp_counts['count'].max()/2
  # ctp_counts['gene'] = ctp_counts['gene'].astype('category').values.reorder_categories(ordered_features)
  ctp_counts['gene'] = ctp_counts['gene'].astype('category').values.reorder_categories(features)
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
        strip_text_y = element_text(size=18, face = 'bold', angle=-90)
      ) +
      scale_x_discrete(labels = lambda list: [re.search(r'^(.*?)_',x).group(1) for x in list]) + 
      coord_flip()
  ).save(img_dir, width=8, height=1+8*len(features), dpi=dpi, 
           limitsize=False, verbose=False)
  return img_dir

def show_features_series_matplotlib(adata, embedding, features=None, pattern=None, sort=True, ascending=True, 
                                          figsize=(6.4,4.8), n_cols=1, dot_size=4, cmap=matplotlib.cm.viridis, **kws):
  
  embedding = embedding.loc[adata.obs_names,]
  
  if not features and pattern:
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()
    
  exp_df = adata[:,features].to_df()
  
  n_rows = int(np.ceil(len(features) / n_cols))
  figsize = (figsize[0]*n_cols, figsize[1]*n_rows)
  fig = plt.figure(figsize=figsize)

  i=1
  for feature in features:
    exp_vec = exp_df[feature]
    if sort:
      exp_vec = exp_vec.sort_values(ascending=ascending)
    embedding_plot = embedding.loc[exp_vec.index]
    
    ax = plt.subplot(n_rows, n_cols, i)
    i = i+1
    plt.scatter(
      x = embedding_plot.iloc[:,0],
      y = embedding_plot.iloc[:,1],
      c = exp_vec,
      cmap = cmap,
      s = dot_size,
      vmin = 0,
      vmax = exp_vec.max()
    )
    plt.xlabel(embedding_plot.columns[0])
    plt.ylabel(embedding_plot.columns[1])
    plt.title(feature, fontsize=16)
    # colorbar
    normalize = matplotlib.colors.Normalize(vmin=0, vmax=exp_vec.max())
    scalarmappable = matplotlib.cm.ScalarMappable(norm=normalize, cmap=cmap)
    scalarmappable.set_array(exp_vec)
    plt.colorbar(scalarmappable, ax=ax)

  if len(features) > 1:
    plt.tight_layout()

  # save fig to a tempory buffer
  buf = BytesIO()
  fig.savefig(buf, format="png")
  # Embed the result in the html output.
  fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
  fig_matplotlib = f'data:image/png;base64,{fig_data}'

  return fig_matplotlib

# In[]: app/widgets/plotFeature

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

num_stages = [ float(re.sub('E','',i)) for i in stage_list]
num_stages = [int(i) if (i % 1 == 0) else i for i in num_stages ]

# In[]: app/widgets/plotSeries

reik_dropdown_featureType_series = html.Div(
    [
        dbc.Label("Feature type"),
        dcc.Dropdown(
            ['Gene'],
            'Gene',
            id="reik_dropdown_featureType_series",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)
reik_input_featureName_series = html.Div(
    [
        dbc.Label("Contains"),
        dbc.InputGroup(
          [
            dmc.Grid(
              children=[
                dmc.Col(dbc.Input(id="reik_input_featureName_series"), span=9),
                dmc.Col(dbc.Button('Plot', id='reik_inputButton_featureName_series_plot', n_clicks=0, color='primary'),
                        span=3),
                dmc.Col(dmc.Text(id='reik_text_seriesGeneNumber_series', color='gray'),
                        span=12),
              ], gutter=3
            )
          ]
        )
    ],
    className="mb-4",
)
reik_textarea_featureLists_series = html.Div(
  [
    dbc.Label("name list:"),
    dbc.Col(
      [
        dbc.Textarea(id = "reik_textarea_featureLists_series",
          placeholder="paste feature names here(seperated by any charactor)\ne.g.  A  B,C,  D\nE##G_H,#I@KK%%G%(U),(V)|(W)\"X\",\"Y\"^Q*I",
          rows=8, className="mb-3",),
        dbc.Button('Plot', id='reik_inputButton_featureLists_series_plot',
                          n_clicks=0, color='primary'),
      ]
    )
  ],
  className="mb-4",
)
reik_dropdown_stage_series = html.Div(
    [
        dbc.Label("Stage"),
        dcc.Dropdown(
            ['E7.5', 'E7.75', 'E8.0', 'E8.5', 'E8.75'],
            'E7.5',
            id="reik_dropdown_stage_series",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

# In[]: app/tabs

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

reik_tab_plotFeatureSeries = dbc.Tab(
  [dbc.Row([
    dbc.Col([
      dbc.Card(
        [
          dbc.Col([
            reik_dropdown_featureType_series,
            reik_dropdown_stage_series,
            reik_input_featureName_series,
            reik_textarea_featureLists_series,
          ])
        ],
        body=True,
        style = {'position': 'fixed', 'width': '30vh'},
        id = 'reik_control_plotFeatureSeries'
      )
    ], width=2),
    dbc.Col([
      dbc.Row(
        [
          dbc.Col(width=2),
          dbc.Col([
            dcc.Graph(figure={}, id="reik_plotFeatureSeries_graph_ctp", style={'height': "40vh", 'width': '80vh'}),
          ],align = "left",className="g-0", width=8),
        ]
      ),
      dbc.Row([
        dbc.Col(
          [
            html.Img(id = 'reik_plotFeatureSeries_img', style = {'width': '160vh'})
          ],
          align = "center", className="g-0", width=12),
      ],),
    ], width=10)
  ])],
  label = "Plot feature(Series)",
  tab_id = "reik_tab_plotFeatureSeries"
)

# In[]: app/all-layout

reik_tabs = dbc.Tabs(
  [reik_tab_plotFeature, reik_tab_plotFeatureSeries],
  active_tab = "reik_tab_plotFeature",  
  id = "reik_tab_plotFeature",
)

layout = dbc.Container(
  [
    html.Div(id='notifications-container-reik'),
    dbc.Row([
      dbc.Col([
        reik_tabs,
      ], width=12)
    ],)
  ],
  fluid=True,
  className="Container-all",
)

# In[]: callbacks/plotFeature

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
  return show_celltype_umap(adata=exp_data[stage], embedding=umap[stage], cmap=celltype_colors, sort=True, orientation='h')

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

# In[]: callbacks/plotSeries

@callback(
  Output('reik_plotFeatureSeries_img', 'src', allow_duplicate=True),
  # Output('reik_plotFeatureSeries_ctpCounts_img', 'src', allow_duplicate=True),
  State('reik_dropdown_featureType_series', 'value'),
  State('reik_input_featureName_series', 'value'),
  Input('reik_inputButton_featureName_series_plot', 'n_clicks'),
  State('reik_dropdown_stage_series', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('reik_inputButton_featureLists_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('reik_inputButton_featureLists_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('reik_inputButton_featureLists_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('reik_textarea_featureLists_series', 'disabled', allow_duplicate=True),True, False),
    (Output('reik_inputButton_featureName_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('reik_inputButton_featureName_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('reik_inputButton_featureName_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('reik_input_featureName_series', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_reik_plotFeature_graphSeries_pattern(featureType, pattern, click, stage):
  if pattern is None:
    raise PreventUpdate

  if click:
    if featureType == 'Gene':
      adata = exp_data[stage]
    elif featureType == 'Regulon':
      adata = auc_mtx[stage]

    img = show_features_series_matplotlib( adata, pattern=pattern, embedding = umap[stage], n_cols=4)
    return img
  else:
    raise PreventUpdate

@callback(
  Output('reik_plotFeatureSeries_img', 'src', allow_duplicate=True),
  # Output('reik_plotFeatureSeries_ctpCounts_img', 'src', allow_duplicate=True),
  Output('notifications-container-reik', 'children'),
  State('reik_dropdown_featureType_series', 'value'),
  State('reik_textarea_featureLists_series', 'value'),
  Input('reik_inputButton_featureLists_series_plot', 'n_clicks'),
  State('reik_dropdown_stage_series', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('reik_inputButton_featureLists_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('reik_inputButton_featureLists_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('reik_inputButton_featureLists_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('reik_textarea_featureLists_series', 'disabled', allow_duplicate=True),True, False),
    (Output('reik_inputButton_featureName_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('reik_inputButton_featureName_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('reik_inputButton_featureName_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('reik_input_featureName_series', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_reik_plotFeature_graphSeries_list(featureType, names, click, stage):

  if names is None:
    raise PreventUpdate

  if click:
    names = re.split(", |,| |\n|\'|\"|#|_|%|$|@|\(|\)|\||^|&", names)
    tmp = [i for i in names if i]
    names = list(set(tmp))

    if featureType == 'Gene':
      adata = exp_data[stage]
    elif featureType == 'Regulon':
      adata = auc_mtx[stage]

    names_out = [i for i in names if (i not in adata.var_names)]
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
    names.sort(key=tmp.index)

    img = show_features_series_matplotlib(adata, features=names, embedding = umap[stage], n_cols=4)
    return img, note
  else:
    raise PreventUpdate

@callback(
  Output('reik_plotFeatureSeries_graph_ctp', 'figure'),
  Input('reik_dropdown_stage_series', 'value'),
)
def update_reik_plotFeatureSeries_ctpGraph(stage):
  fig = show_celltype_umap_series(exp_data[stage], embedding = umap[stage], cmap=celltype_colors, sort=True)
  return fig

@callback( # update series-gene's number
  Output('reik_text_seriesGeneNumber_series', 'children'),
  Input('reik_dropdown_featureType_series', 'value'),
  Input('reik_dropdown_stage_series', 'value'),
  Input('reik_input_featureName_series', 'value'),
)
def update_spatial_text_seriesGeneNumber_series(featureType, stage, pattern):
  if featureType == 'Gene':
    adata = exp_data[stage]
  elif featureType == 'Regulon':
    adata = auc_data[stage]
  
  features = [i for i in adata.var_names if re.match(pattern, i)]
  str = f'{len(features)} genes'
  return str

