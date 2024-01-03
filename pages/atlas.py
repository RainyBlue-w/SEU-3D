#!/usr/bin/env python
# coding: utf-8

# In[1]:

import dash
dash.register_page(__name__)


from dash import Dash, dcc, html, dash_table, Input, Output, callback, no_update, State, Patch, DiskcacheManager
from dash.exceptions import PreventUpdate
import plotly.express as px
import plotly
import dash_bootstrap_components as dbc
import scanpy as sc
import os
import pandas as pd
import numpy as np
import loompy as lp

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

data_dir = "/rad/wuc/dash_data/atlas/"

# color palette
celltype = ["ExE endoderm", "Epiblast", "ExE ectoderm", "Parietal endoderm", 
  "Intermediate mesoderm", "Gut", "Caudal epiblast", "Haematoendothelial progenitors", 
  "Visceral endoderm", "PGC", "Somitic mesoderm", "Def. endoderm", "Mesenchyme", 
  "Blood progenitors 2", "Nascent mesoderm", "Paraxial mesoderm", "Notochord", 
  "Mixed mesoderm", "Blood progenitors 1", "Rostral neurectoderm", 
  "Anterior Primitive Streak", "Surface ectoderm", "Primitive Streak", "ExE mesoderm", 
  "Forebrain Midbrain Hindbrain", "Pharyngeal mesoderm", "Erythroid1", "Cardiomyocytes",
   "Caudal Mesoderm", "Endothelium", "Caudal neurectoderm", "Spinal cord", "NMP",
    "Erythroid3", "Erythroid2", "Allantois", "Neural crest"]
celltype_color = ["#7f6874", "#635547", "#989898", "#1a1a1a", "#139992", "#ef5a9d",
  "#9e6762", "#fbbe92", "#f6bfcb", "#facb12", "#005579", "#f397c0", "#cc7818", "#c9a997",
  "#c594bf", "#8db5ce", "#0f4a9c", "#dfcde4", "#f9decf", "#65a83e", "#c19f70", "#f7f79e",
  "#dabe99", "#8870ad", "#647a4f", "#c9ebfb", "#c72228", "#b51d8d", "#3f84aa", "#ff891c",
  "#354e23", "#cde088", "#8ec792", "#ef4e22", "#f79083", "#532c8a", "#c3c388"]
timepoint = ["E6.5", "E6.75" ,"E7.0", "E7.25", "E7.5", "E7.75", "E8.0",
  "E8.25", "E8.5", "mixed_gastrulation", "E8.75"]
timepoint_color = ["#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#ffffbf", 
  "#e6f598", "#abdda4", "#66c2a5", "#3288bd", "#a9a9a9", "#B53588"]

color_palette = dict(zip(celltype, celltype_color))
color_palette.update(dict(zip(timepoint, timepoint_color)))


# In[2]:


import os
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import celloracle as co
co.__version__


# Code below shows how the anndata with diffgene results was generated:

# In[3]:


# adata = sc.read_h5ad(data_dir+ 'Time_series_tSNE/atlas.var.h5ad')
# celltype=adata.obs["celltype"].unique()
# for i in range(len(celltype)):
    # sc.tl.rank_genes_groups(adata, 'celltype', groups=celltype.delete(i).to_numpy(), reference=celltype[i], method="wilcoxon", key_added="diff_ref:"+celltype[i])
# adata.write_h5ad(data_dir+ 'Time_series_tSNE/atlas.var.diff_res.h5ad')


# In[7]:
def get_pdf(stage):
    pdf = pd.read_excel((data_dir+'Time_series_tSNE/tSNE_excel/%s.xlsx' % (stage)), index_col=0)
    metadata = pd.read_csv(data_dir+'Time_series_tSNE/metadata/meta_'+stage+'.tab', delimiter='\t').loc[pdf.index]
    return pdf


# prepare data
stage_list = ['E7.0',
              'E7.25',
              'E7.5',
              'E7.75',
              'E8.0',
              'E8.25',
              'E8.5',
             ]

pdf = dict()
df_tSNE = None
for i in stage_list:
    pdf[i] = get_pdf(i)
    if df_tSNE is None:
        df_tSNE = pdf[i]
    else:
        df_tSNE = pd.concat([df_tSNE, pdf[i]], axis=0)

adata = sc.read_h5ad(data_dir+ 'Time_series_tSNE/atlas.var.h5ad')
all(df_tSNE.index.isin(adata.obs.index))
adata = adata[df_tSNE.index,]
adata.obs = pd.concat([adata.obs, df_tSNE.loc[adata.obs.index,]], axis=1)


# In[8]:


lf = lp.connect('/rad/wuc/dash_data/atlas/pyscenic_output.loom', mode='r', validate=False )
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

def show_gene(gene, adata, embedding=None, cmap = None):
    if cmap is None:
        cmap = [(0.00, "rgb(150,150,150)"),
                (50/256, "rgb(217, 217, 217)"),
                (50/256, "rgb(191, 217, 237)"),
                (1.00, "rgb(8, 48, 107)")
                ]
    if embedding is None:
        embedding = adata.obs[['tSNE-1', 'tSNE-2']]
    pdf = pd.DataFrame(np.array(embedding), 
                       index=adata.obs_names, 
                       columns=['tSNE-1', 'tSNE-2'])

    pdf[gene] = adata.to_df()[gene].tolist()
    plot = px.scatter(
    		data_frame = pdf,
    		x = 'tSNE-1', y = 'tSNE-2', color = gene,
    		color_continuous_scale = cmap
    	)
    plot.update_layout(showlegend=False)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=4.5,
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

def show_regulon(regulon, auc_mtx, stage, cmap = None):
    if cmap is None:
        cmap = [(0.00, "rgb(150,150,150)"),
                (50/256, "rgb(217, 217, 217)"),
                (50/256, "rgb(191, 217, 237)"),
                (1.00, "rgb(8, 48, 107)")
                ]
    tmp = auc_mtx.index.isin(adata[adata.obs.stage==stage].obs_names)
    cell_id = auc_mtx.loc[tmp,:].index
    auc = auc_mtx.loc[cell_id,:]
    embedding = adata[cell_id].obs[['tSNE-1', 'tSNE-2']]
    pdf = pd.DataFrame(np.array(embedding),
                      index = cell_id,
                      columns=['tSNE-1', 'tSNE-2'])
    pdf[regulon] = auc[regulon].tolist()
    plot = px.scatter(
        data_frame = pdf,
        x = 'tSNE-1', y = 'tSNE-2', color = regulon,
        color_continuous_scale = cmap
    )
    plot.update_layout(showlegend=False)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=4.5,
                    marker_opacity=1)
    plot.update_layout(
        margin=dict(l=10, r=10, t=10, b=10),
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
            color_discrete_map = cmap,
    )
    # plot.update_layout(showlegend=False)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
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

def get_regulon_targets(regulon):
    return pd.DataFrame(regulons[regulon])

def get_celltype(pdf):
    return pdf["celltype"].unique()


# In[10]:


# define components


dropdown_gene = html.Div(
    [
        dbc.Label("Gene to display"),
        dcc.Dropdown(
            adata.var_names,
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
        dbc.Label("Regulon's target gene"),
        dash_table.DataTable(sort_action="native", page_action='native',
                             page_current= 0, page_size= 10,
                             id = "target_gene",fill_width=True, style_table={'overflowX': 'auto'},
                                style_cell={
                                'padding-right': '30px',
                                'padding-left': '10px',
                                'text-align': 'center',
                                'marginLeft': 'auto',
                                'marginRight': 'auto'
                                },
                             filter_action="native", 
                            # filter_options={'case':'insensitive'}, 
                            ),
    ],
    className="mb-4",
)

dropdown_regulon = html.Div(
    [
        dbc.Label("Regulon to display"),
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

dropdown_diffgene_ctp1 = html.Div(
    [
        dbc.Label("celltype 1"),
        dcc.Dropdown(
            [],
            id="diffgene_ctp1",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)

dropdown_diffgene_ctp2 = html.Div(
    [
        dbc.Label("celltype 2"),
        dcc.Dropdown(
            [],
            id="diffgene_ctp2",
            clearable=False,
            searchable=True,
        ),
    ],
    className="mb-4",
)


dropdown_markergene_ctp = html.Div(
    [
        dbc.Label("celltype"),
        dcc.Dropdown(
            [],
            id="markergene_ctp",
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

stage_diffgene = dcc.Slider(0, 6,
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
    id="stage_diffgene"
)

stage_markergene = dcc.Slider(0, 6,
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
    id="stage_markergene"
)

tsne_controls = dbc.Card(
    [
        # dropdown_projType,
        dropdown_gene,
        dropdown_regulon,
        target_gene,
    ],
    body=True,
)

diffgene_controls = dbc.Card(
    [
        dropdown_diffgene_ctp1,
        dropdown_diffgene_ctp2,
    ],
    body=True,
)

markergene_controls = dbc.Card(
    [
        dropdown_markergene_ctp,
    ],
    body=True,
)

# controls layout
controls = dbc.Col(
    tsne_controls,
    id = "controls",
    style = {'position': 'fixed'},
    width=2
)

# tabs layout
tab_tsne = dbc.Tab(
    [
        dbc.Col(
            [
                dbc.Row(
                    [
                        dbc.Col(
                          dbc.Row(
                            [
                              # html.H5('Celltype',style={'textAlign': 'center'}),
                              dcc.Graph(figure={}, id="tsne_ctp",style={'height': "75vh"})
                            ], 
                            className="g-0",
                          ),
                          width=6, xxl=4,
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
                                width=12, xxl=6,
                              ),
                              dbc.Col(
                                [
                                  html.H5("Regulon", style={'textAlign': 'center'}, id='regulon_title'),
                                  dcc.Graph(figure={}, id="tsne_regulon", style={'height': "60vh"}),
                                ],
                                width=12, xxl=6,
                              )
                              ],
                              className="g-0",
                            )
                          ],
                          width=6, xxl=8,
                        ),
                    ],
                    align = "center",
                    className="g-0",
                ),
                stage_tsne,
            ], 
        )  
    ], 
    label="tsne",
    tab_id = "tab_tsne"
)

df = pd.read_csv('/home/wuc/dashapps/multi_page/resources/solar.csv')
tab_diffgene = dbc.Tab(
    [
        dbc.Row(
            [
                dbc.Col(
                    dcc.Graph(figure={}, id="diffgene_tsne"),
                    width=6
                ),
                dbc.Col(
                    dash_table.DataTable(df.to_dict('records'), [{"name": i, "id": i} for i in df.columns],
                                        sort_action="native", sort_mode="multi", row_selectable="single"),
                    width=6
                ),
            ]
        ),
        stage_diffgene,
    ],
    label = "diff gene",
    tab_id = "tab_diffgene"
)

tab_markergene = dbc.Tab(
    [
        dbc.Row(
            [
                dbc.Col(
                    dcc.Graph(figure={}, id="markergene_tsne"),
                    width=6,
                ),
                dbc.Col(
                    dash_table.DataTable(df.to_dict('records'), [{"name": i, "id": i} for i in df.columns],
                                        sort_action="native", sort_mode="multi", row_selectable="single"),
                    width=6
                ),
            ],
        ),
        stage_markergene,
    ],
    label = "marker gene",
    tab_id = "tab_markergene"
)

tabs = dbc.Card(
    dbc.Tabs(
        [tab_tsne, tab_diffgene, tab_markergene],
        active_tab = "tab_tsne",  
        id = "tabs",
    ),
)

# all layout
layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [controls],
                    width=2,
                ),
                dbc.Col(
                    [tabs],
                   width=10
                )
            ],
        ),
    ],
    fluid=True,
    className="dbc",
)

@callback(
    Output("controls", "children"),
    Input("tabs", "active_tab"), 
)
def update_controls(tab):
    if tab == "tab_tsne":
        return tsne_controls
    elif tab == "tab_diffgene":
        return diffgene_controls
    elif tab == "tab_markergene":
        return markergene_controls

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
    plot_gene = show_gene(gene, adata[adata.obs.stage==stage], 
                          cmap = 'Inferno')
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
    plot_regulon = show_regulon(regulon, auc_mtx, stage, 
                               cmap = 'Inferno')
    return regulon, plot_regulon

@callback(
  Output("gene", 'value'),
  Input("target_gene", 'active_cell'),
  State("regulon", "value"),
)
def update_tsne_byTargetGene(active_cell, regulon):
  # if not active_cell:
  #   raise PreventUpdate
  row = active_cell['row_id']
  df = get_regulon_targets(regulon)
  return df.iloc[row,0]

@callback(
    Output("diffgene_ctp1", "options"),
    Output("diffgene_ctp1", "value"),
    Output("diffgene_ctp2", "options"),
    Output("diffgene_ctp2", "value"),
    Input("stage_diffgene", "value")
)
def update_diffgene_ctp(value):
    stage = stage_list[value]
    ctp = pdf[stage]["celltype"].unique()
    return ctp, ctp[1], ctp, ctp[1]

@callback(
    Output("markergene_ctp", "options"),
    Output("markergene_ctp", "value"),
    Input("stage_markergene", "value")
)
def update_diffgene_ctp(value):
    stage = stage_list[value]
    ctp = pdf[stage]["celltype"].unique()
    return ctp, ctp[1]


# In[11]:





# In[ ]:




