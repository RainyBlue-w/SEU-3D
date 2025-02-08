#!/usr/bin/env python
# coding: utf-8

import dash
dash.register_page(__name__)

# In[1]:

from dash import Dash, dcc, html, dash_table, Input, Output, callback, no_update, State, Patch, DiskcacheManager
from dash.exceptions import PreventUpdate
import plotly.express as px
import dash_bootstrap_components as dbc
import scanpy as sc
import os
import pandas as pd
import numpy as np
import loompy as lp
import os
import numpy as np
import pandas as pd
import scanpy as sc
from dash_iconify import DashIconify
from concurrent import futures
from PIL import Image
import dash_mantine_components as dmc
from plotnine import *
import re
import diskcache
import matplotlib 
import matplotlib.pyplot as plt
import base64
from io import BytesIO
import matplotlib.colors as mcolors
from utils.ctp_colormap import ctp_cmap
matplotlib.use('agg')

data_dir = "/data1/share/omics-viewer/atlas-CP/"

background_callback_manager = DiskcacheManager(diskcache.Cache(data_dir+"cache"))


# color palette
color_palette = ctp_cmap

timepoint = ["E6.5", "E6.75" ,"E7.0", "E7.25", "E7.5", "E7.75", "E8.0",
  "E8.25", "E8.5", "mixed_gastrulation", "E8.75"]
timepoint_color = ["#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", 
  "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#a9a9a9", "#B53588"]

color_palette.update(dict(zip(timepoint, timepoint_color)))

from matplotlib.colors import LinearSegmentedColormap
color_exp = [
        (0.00, "#BFBFBF"),
        (0.75, "#225EA8"),
        (1.00, "#000000")
    ]
cmap_exp = LinearSegmentedColormap.from_list('custom_exp', color_exp)

color_auc = [(0.00, "#BFBFBF"),
              (0.75, "#d53e4f"),
              (1.00, "#000000")
            ]
cmap_auc = LinearSegmentedColormap.from_list('custom_auc', color_auc)
        

# In[7]:
def get_pdf(stage):
    pdf = pd.read_excel((data_dir+'Time_series_tSNE/tSNE_excel/%s.xlsx' % (stage)), index_col=0)
    return pdf

# prepare data
stage_list = ['E7.0',
              'E7.25',
              'E7.5',
              'E7.75',
              'E8.0',
              'E8.25',
              'E8.5']

pdf = dict()
df_tSNE = None
for i in stage_list:
    pdf[i] = get_pdf(i)
    if df_tSNE is None:
        df_tSNE = pdf[i]
    else:
        df_tSNE = pd.concat([df_tSNE, pdf[i]], axis=0)

exp_data = sc.read_h5ad(data_dir+ 'backedModeData/atlasCP.h5ad', backed='r')

# exp_data = sc.read_h5ad(data_dir+ 'Time_series_tSNE/atlas.var.h5ad', backed='r')
# exp_data = exp_data.raw.to_adata()
# exp_data = exp_data[df_tSNE.index,]
# exp_data.obs[['tSNE-1', 'tSNE-2']] = df_tSNE.loc[exp_data.obs.index][['tSNE-1', 'tSNE-2']]
# exp_data.obs.index.name = None
# In[8]:

lf = lp.connect(data_dir+'pyscenic_output.loom', mode='r', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).items():
    regulons[i] =  list(r[r==1].index.values)
lf.close()


# In[9]:

# functions
def convert_color(s):
    """
    Convert a string hexadecimal representation of a color into a tuple
    RBG where each element of the tuple is between 0.0 to 1.0

    :param s: (str)
    :return: (tuple) With 3 floats representing the color in RGB

    :rtype : tuple

    Examples:
    >>> convert_color('FF5500')
    (1.0, 0.3333333333333333, 0.0)

    """
    return float(int(s[:2], 16)) / 255, float(int(s[2:4], 16)) / 255, float(int(s[4:6], 16)) / 255 

def show_gene(gene, adata, embedding=None, cmap = color_exp):
    if embedding is None:
        embedding = adata.obs[['tSNE-1', 'tSNE-2']]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['tSNE-1', 'tSNE-2'])

    pdf[gene] = adata.to_df()[gene].tolist()
    pdf = pdf.sort_values(by=gene, ascending=True)
    plot = px.scatter(
      data_frame = pdf,
      x = 'tSNE-1', y = 'tSNE-2', color = gene,
      color_continuous_scale = cmap
    )
    plot.update_layout(showlegend=False)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=5,
                    marker_opacity=1)
    plot.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        plot_bgcolor = '#ffffff',
        legend = dict(
            title='Expression'
        ),
        uirevision='constant',
        legend_itemsizing = 'constant'
    )
    return plot

def show_regulon(regulon, auc_mtx, stage, cmap = color_auc):
    tmp = auc_mtx.index.isin(exp_data[exp_data.obs.stage==stage].obs_names)
    cell_id = auc_mtx.loc[tmp,:].index
    auc = auc_mtx.loc[cell_id,:]
    embedding = exp_data[cell_id].obs[['tSNE-1', 'tSNE-2']]
    pdf = pd.DataFrame(np.array(embedding),
                      index = cell_id,
                      columns=['tSNE-1', 'tSNE-2'])
    pdf[regulon] = auc[regulon].tolist()
    pdf = pdf.sort_values(by=regulon, ascending=True)
    plot = px.scatter(
        data_frame = pdf,
        x = 'tSNE-1', y = 'tSNE-2', color = regulon,
        color_continuous_scale = cmap
    )
    plot.update_layout(showlegend=False)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=5,
                    marker_opacity=1)
    plot.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        plot_bgcolor = '#ffffff',
        legend = dict(
            title='AUCell Score'
        ),
        uirevision='constant',
        legend_itemsizing = 'constant'
    )
    return plot

def show_celltype(pdf, cmap):
    plot = px.scatter(
            data_frame = pdf,
            x='tSNE-1', y='tSNE-2', color='celltype',
            color_discrete_map = cmap
    )
    # plot.update_layout(showlegend=False)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=5,
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
        )
    )
    return plot

def show_celltype_series(pdf, cmap, figwidth=800, figheight=300):
    plot = px.scatter(
            data_frame = pdf,
            x='tSNE-1', y='tSNE-2', color='celltype',
            color_discrete_map = cmap,
            width=figwidth, height=figheight,
    )
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=3, marker_opacity=1)
    plot.update_layout(
      margin=dict(l=0, t=0, b=0),
      plot_bgcolor = '#ffffff',
      title = '',
      legend=dict(
              title='',
              orientation='h',
              yanchor='middle',
              xanchor='right',
              y=0.5,
              x=2,
      ),
      uirevision='constant',
      legend_itemsizing = 'constant'
    )
    return plot

def get_regulon_targets(regulon):
    return pd.DataFrame(regulons[regulon])

def get_celltype(pdf):
    return pdf["celltype"].unique()

def show_features_atlas_regularExp(adata, stage, odir, featureType, embedding, pattern=None, features=None, sort=False, ascending=True, dpi=100, **kws):
  embedding = embedding.loc[adata.obs_names,:]

  if featureType == 'Regulon':
    colors = color_auc
  else: 
    colors = color_exp 
  
  if not features:
    
    img_dir = data_dir+'tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()
    # features = [i for i in features if i in genes_min_pval.loc[stage].index.to_list()]

  else:
    img_dir = data_dir+'tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)

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
    scale_color_gradientn(colors = [i[1] for i in colors], values = [i[0] for i in colors]) +
    theme(
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

def show_featuresCtpcounts_atlas_regularExp(adata, stage, odir, featureType, embedding, pattern=None, features=None, cmap=color_palette, sort=False, ascending=True, dpi=150, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  if not features:
    img_dir = data_dir+'tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()

  else:
    img_dir = data_dir+'tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)

  ctp_counts = {}
  for gene in features:
    df = adata[:,gene].to_df()
    df = df[df[gene] > 0]
    counts = pd.DataFrame(adata.obs['celltype'][df.index].value_counts())
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
        strip_text_y = element_text(size=16, face = 'bold', angle=-90)
      ) +
      scale_x_discrete(labels = lambda list: [re.search(r'^(.*?)_',x).group(1) for x in list]) + 
      coord_flip()
  ).save(img_dir, width=8, height=1+8*len(features), dpi=dpi, 
           limitsize=False, verbose=False)
  return img_dir

def show_features_series_matplotlib_atlas(adata, stage, embedding=None, features=None, pattern=None, sort=True, ascending=True, 
                                          figsize=(6.4,4.8), n_cols=1, dot_size=4, cmap=cmap_exp, **kws):
  
  if not features and pattern:
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()

  adata = adata[adata.obs.stage==stage, features]

  embedding = embedding.loc[adata.obs_names,]

  exp_df = adata.to_df()
  
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


# In[10]:

dropdown_gene = html.Div(
    [
        dbc.Label("Gene to Display"),
        dcc.Dropdown(
            exp_data.var_names,
            'Gata6',
            id="gene",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

target_gene = html.Div(
    [
        dbc.Label("Regulon's Target Gene"),
        dash_table.DataTable(sort_action="native", page_action='native',
                             page_current= 0, page_size= 9,
                             id = "target_gene",fill_width=True, style_table={'overflowX': 'auto'},
                                style_cell={
                                'text-align': 'center',
                                },
                             filter_action="native", 
                            # filter_options={'case':'insensitive'}, 
                            ),
    ],
    className="mb-4",
)

dropdown_regulon = html.Div(
    [
        dbc.Label("Regulon to Display"),
        dcc.Dropdown(
            auc_mtx.columns,
            'Gata6(+)',
            id="regulon",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

stage_tsne = dcc.Slider(0, 6,
    step = None,
    marks={
        0: 'E7.0',
        1: 'E7.25',
        2: 'E7.5',
        3: 'E7.75',
        4: 'E8.0',
        5: 'E8.25',
        6: 'E8.5',
    },
    value=5,
    id="stage_tsne"
)

atlas_dropdown_featureType_series = html.Div(
    [
        dbc.Label("Feature Type"),
        dcc.Dropdown(
            ['Gene'],
            'Gene',
            id="atlas_dropdown_featureType_series",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)
atlas_input_featureName_series = html.Div(
    [
        dbc.Label("Contains"),
        dbc.InputGroup(
          [
            dmc.Grid(
              children=[
                dmc.GridCol(dbc.Input(id="atlas_input_featureName_series"), span=9),
                dmc.GridCol(dbc.Button('Plot', id='atlas_inputButton_featureName_series_plot', 
                                   n_clicks=0, color='primary'),
                        span=3),
                dmc.GridCol(dmc.Text(id='atlas_text_seriesGeneNumber_series', c='gray'),
                        span=12),
              ], gutter=3
            )
          ]
        )
    ],
    className="mb-4",
)
atlas_textarea_featureLists_series = html.Div(
  [
    html.Div(
       [
          dbc.Label("Name List"),
          dbc.Button('Plot', id='atlas_inputButton_featureLists_series_plot',
                          n_clicks=0, color='primary', style={'float': 'right', 'margin-bottom': '1vh'}),
       ], className="mb-1"
    ),
    dbc.Col(
      [
        dbc.Textarea(id = "atlas_textarea_featureLists_series",
          placeholder="paste feature names here(seperated by any charactor)\ne.g.  A  B,C,  D\nE##G_H,#I@KK%%G%(U),(V)|(W)\"X\",\"Y\"^Q*I",
          rows=8, className="mb-3",)
      ]
    )
  ],
  className="mb-4",
)
atlas_dropdown_stage_series = html.Div(
    [
        dbc.Label("Stage"),
        dcc.Dropdown(
            ['E7.0', 'E7.25', 'E7.5', 'E7.75', 'E8.0', 'E8.25', 'E8.5'],
            'E7.0',
            id="atlas_dropdown_stage_series",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

series_controls = dbc.Card(
    dbc.Col(
        [
            atlas_dropdown_featureType_series,
            atlas_dropdown_stage_series,
            atlas_input_featureName_series,
            atlas_textarea_featureLists_series,
        ]
    ),
    body=True,
    style = {'width': '16.67vw', 'height':'83vh'},
    id = 'atlas_control_plotFeatureSeries'
)

tsne_controls = dbc.Card(
    [
        # dropdown_projType,
        
    ],
    body=True,
)

# controls layout
tsne_controls = dbc.Card(
   dbc.Col([
        dropdown_gene,
        dropdown_regulon,
        target_gene,
   ]),
    id = "tsne_controls",
    body=True,
    style = {'position': 'fixed', 'width': '16.67vw', 'height':'100vh'},
)

# tabs layout
tab_tsne = dbc.Tab(
    [
        dbc.Row([
           dbc.Col(
            tsne_controls,
            width = 2
           ),
           dbc.Col(
            [
                dbc.Row(
                    [
                        dbc.Col(
                          dbc.Row(
                            [
                              # html.H5('Celltype',style={'textAlign': 'center'}),
                              dcc.Graph(figure={}, id="tsne_ctp",style={'height': "76vh"})
                            ], 
                          ),
                          width=4
                        ),
                        dbc.Col(
                          [
                            dbc.Row(
                              [
                                dbc.Col(
                                  [
                                    html.H5("Gene", style={'textAlign': 'center'}, id='gene_title'),
                                    dcc.Graph(figure={}, id="tsne_gene", style={'height': "60vh"}),
                                  ],
                                  width=6
                                ),
                                dbc.Col(
                                  [
                                    html.H5("Regulon", style={'textAlign': 'center'}, id='regulon_title'),
                                    dcc.Graph(figure={}, id="tsne_regulon", style={'height': "60vh"}),
                                  ],
                                  width=6
                                )
                              ],
                            )
                          ]
                        ),
                    ],
                    align = "center",
                    className="g-0",
                ),
                stage_tsne,
            ], width = 10
          )
        ])], 
    label="tsne",
    tab_id = "tab_tsne"
)

tab_series = dbc.Tab(
  [
    dbc.Row(
       [
          dbc.Col([
            series_controls
          ], width=2),
          dbc.Col([
            dbc.Row([
              dbc.Col(width=2),
              dbc.Col([
                dcc.Graph(figure={}, id="atlas_plotFeatureSeries_graph_ctp", style={'height': "40vh", 'width': '80vh'}),
              ],align = "center",className="g-0", width=8),
            ]),
            dbc.Row([
              dbc.Col(
                [
                  html.Img(id = 'atlas_plotFeatureSeries_img', style = {'width': '80vw', 'margin-top': '5vh', 'margin-left':'3vw'})
                ],
                align = "center", className="g-0", width=12),
            ],),
          ], width=10)
       ]
    )
  ],
  label = "Plot Feature(Series)",
  tab_id = "tab_series"
)

tabs = dbc.Card(
    dbc.Tabs(
        [tab_tsne, tab_series],
        active_tab = "tab_tsne",  
        id = "tabs",
    ), 
    style={'height':'90vh'}
)

# all layout
layout = dbc.Container(
    [
        html.Div(id='notifications-container-atlas'),
        tabs
    ],
    fluid=True,
    className="dbc",
)

@callback(
    Output("target_gene", "data"),
    Output("target_gene", "columns"),
    Input("regulon", "value")
)
def update_target_gene(regulon):
    df = get_regulon_targets(regulon)
    df.columns = ['Genes']
    df = df.reset_index().rename(columns={"index": "id"})
    return df.to_dict('records'), [{"name": i, "id": i, "deletable": True} for i in df.columns if i != 'id']

@callback(
    Output('gene_title', 'children'),
    Output("tsne_gene", "figure"),

    Input("stage_tsne", "value"),
    Input("gene", "value"),
)
def update_tsne_gene(value, gene):
    stage = stage_list[value]
    plot_gene = show_gene(gene, exp_data[exp_data.obs.stage==stage])
    return gene, plot_gene

@callback(
    Output("tsne_ctp", "figure"),
    Input("stage_tsne", "value"),
)
def update_tsne_ctp(value):
    stage = stage_list[value]
    plot_ctp = show_celltype(pdf[stage], color_palette)
    return plot_ctp

@callback(
    Output('regulon_title', 'children'),
    Output("tsne_regulon", "figure"),
    
    Input("stage_tsne", "value"),
    Input("regulon", "value"),
)
def update_tsne_regulon(value, regulon):
    stage = stage_list[value]
    plot_regulon = show_regulon(regulon, auc_mtx, stage)
    return regulon, plot_regulon

@callback(
  Output("gene", 'value'),
  Input("target_gene", 'active_cell'),
  State("regulon", "value"),
)
def update_tsne_byTargetGene(active_cell, regulon):
  if not active_cell:
    raise PreventUpdate
  row = active_cell['row_id']
  df = get_regulon_targets(regulon)
  return df.iloc[row,0]

@callback(
  Output('atlas_plotFeatureSeries_img', 'src', allow_duplicate=True),
  # Output('atlas_plotFeatureSeries_ctpCounts_img', 'src', allow_duplicate=True),
  State('atlas_dropdown_featureType_series', 'value'),
  State('atlas_input_featureName_series', 'value'),
  Input('atlas_inputButton_featureName_series_plot', 'n_clicks'),
  State('atlas_dropdown_stage_series', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('atlas_inputButton_featureLists_series_plot', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('atlas_inputButton_featureLists_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('atlas_inputButton_featureLists_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlas_textarea_featureLists_series', 'disabled', allow_duplicate=True),True, False),
    (Output('atlas_inputButton_featureName_series_plot', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('atlas_inputButton_featureName_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('atlas_inputButton_featureName_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlas_input_featureName_series', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_atlas_plotFeature_graphSeries_pattern(featureType, pattern, click, stage):
  if pattern is None:
    raise PreventUpdate

  if click:
    if featureType == 'Gene':
      embedding = exp_data[exp_data.obs.stage==stage].obs[['tSNE-1', 'tSNE-2']]
      adata = exp_data

    # elif featureType == 'Regulon':
    #   adata = auc_mtx[stage]
    else:
        raise PreventUpdate

    img = show_features_series_matplotlib_atlas(adata, stage, embedding=embedding, n_cols=4, pattern=pattern)
    return img
  else:
    raise PreventUpdate

@callback(
  Output('atlas_plotFeatureSeries_img', 'src', allow_duplicate=True),
  Output('notifications-container-atlas', 'children'),
  State('atlas_dropdown_featureType_series', 'value'),
  State('atlas_textarea_featureLists_series', 'value'),
  Input('atlas_inputButton_featureLists_series_plot', 'n_clicks'),
  State('atlas_dropdown_stage_series', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('atlas_inputButton_featureLists_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('atlas_inputButton_featureLists_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('atlas_inputButton_featureLists_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlas_textarea_featureLists_series', 'disabled', allow_duplicate=True),True, False),
    (Output('atlas_inputButton_featureName_series_plot', 'children', allow_duplicate=True), 'Loading', 'Plot'),
    (Output('atlas_inputButton_featureName_series_plot', 'disabled', allow_duplicate=True), True, False),
    (Output('atlas_inputButton_featureName_series_plot', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlas_input_featureName_series', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_atlas_plotFeature_graphSeries_list(featureType, names, click, stage):

  if names is None:
    raise PreventUpdate

  if click:
    names = re.split(", |,| |\n|\'|\"|#|_|%|$|@|\(|\)|\||^|&", names)
    tmp = [i for i in names if i]
    names = list(set(tmp))

    if featureType == 'Gene':
      embedding = exp_data[exp_data.obs.stage==stage].obs[['tSNE-1', 'tSNE-2']]
      adata = exp_data
    # elif featureType == 'Regulon':
    #   adata = auc_mtx[stage]

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

    img = show_features_series_matplotlib_atlas(adata, stage, embedding=embedding, n_cols=4, features=names)
    return img, note
  else:
    raise PreventUpdate

@callback(
  Output('atlas_plotFeatureSeries_graph_ctp', 'figure'),
  Input('atlas_dropdown_stage_series', 'value'),
)
def update_atlas_plotFeatureSeries_ctpGraph(stage):
  fig = show_celltype_series(pdf[stage], color_palette)
  return fig

@callback( # update series-gene's number
  Output('atlas_text_seriesGeneNumber_series', 'children'),
  Input('atlas_dropdown_featureType_series', 'value'),
  Input('atlas_dropdown_stage_series', 'value'),
  Input('atlas_input_featureName_series', 'value'),
)
def update_spatial_text_seriesGeneNumber_series(featureType, stage, pattern):
  if not pattern:
     return ''
  if featureType == 'Gene':
    adata = exp_data[exp_data.obs.stage==stage]
#   elif featureType == 'Regulon':
#     adata = auc_data[stage]
  
  features = [i for i in adata.var_names if re.match(pattern, i)]
  str = f'{len(features)} genes'
  return str