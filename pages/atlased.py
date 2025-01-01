#!/usr/bin/env python
# coding: utf-8

import dash
dash.register_page(__name__)

# In[1]:

from dash import dcc, html, DiskcacheManager
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from plotnine import *
import diskcache
import matplotlib
matplotlib.use('agg')
from callbacks.atlased import *

background_callback_manager = DiskcacheManager(diskcache.Cache('/data1/share/omics-viewer/atlas-Endoderm/cache'))

endoderm_gene_dropdown = html.Div(
  [
    dbc.Label("Gene to Display"),
    dcc.Dropdown(
      geneList,
      geneList[3],
      id="endoderm_gene_dropdown",
      clearable=False,
      searchable=True,
    ),
  ], className="mb-4",
)

endoderm_markerSize_dropdown = html.Div(
  [
    dbc.Label("Marker Size"),
    dcc.Dropdown(
      [1,2,3,4,5,6,7,8,9,10],
      5,
      id="endoderm_markerSize_dropdown",
      clearable=False,
      searchable=True,
    ),
  ], className="mb-4",
)

endoderm_stage_slider = html.Div(
  dcc.Slider(
    0, len(stageDict)-1,
    step = None,
    marks=stageDict,
    value=2,
    id="endoderm_stage_slider"
  ), style={'width': '97%'}
)

endoderm_featureType_series_dropdown = html.Div(
  [
    dbc.Label("Feature Type"),
    dcc.Dropdown(
      ['celltype', 'stage', 'endo_gutCluster'],
      'celltype',
      id="endoderm_featureType_series_dropdown",
      clearable=False,
      searchable=True,
    ),
  ], className="mb-4",
)

endoderm_stage_series_dropdown = html.Div(
  [
    dbc.Label("Stage"),
    dcc.Dropdown(
      stageValues,
      stageValues[3],
      id="endoderm_stage_series_dropdown",
      clearable=False,
      searchable=True,
    ),
  ], className="mb-4",
)

endoderm_marker_gene_dropdown = html.Div(
  [
    dbc.Label("Marker Gene"),
    dmc.Grid(
      [
          dmc.GridCol(dcc.Dropdown(
            id="endoderm_marker_gene_dropdown",
            clearable=False,
            searchable=True,
          ), span=9),
          dmc.GridCol(dbc.Button('Get', id='endoderm_marker_gene_button', n_clicks=0, color='primary'), span=3),
      ], gutter=3
    )
  ], className="mb-4",
)

endoderm_featureName_series_input = html.Div(
  [
    dbc.Label("Contains"),
    dbc.InputGroup(
      [
        dmc.Grid(
          children=[
            dmc.GridCol(dbc.Input(id="endoderm_featureName_series_input"), span=9),
            dmc.GridCol(dbc.Button('Plot', id='endoderm_featureName_series_inputButton', n_clicks=0, color='primary'), span=3),
            dmc.GridCol(dmc.Text(id='endoderm_geneNumber_series_text', c='gray'), span=12),
          ], gutter=3
        )
      ]
    )
  ], className="mb-4",
)

endoderm_featureList_series_textarea = html.Div(
  [
    html.Div(
      dmc.Grid(
        [
          dmc.GridCol(dbc.Label("Name List"), span=9),
          dmc.GridCol(dbc.Button('Plot', id='endoderm_featureList_series_inputButton', n_clicks=0, color='primary'), span=3)
        ], gutter=3
      ), className="mb-2"
    ),
    html.Div(
      dbc.Textarea(id = "endoderm_featureList_series_textarea",
        placeholder="paste feature names here(seperated by any charactor)\ne.g.  A  B,C,  D\nE##G_H,#I@KK%%G%(U),(V)|(W)\"X\",\"Y\"^Q*I",
        rows=8, className="mb-3",)
    )
  ], className="mb-4",
)

endoderm_series_controls = dbc.Card(
  dbc.Col(
    [
      endoderm_featureType_series_dropdown,
      endoderm_stage_series_dropdown,
      endoderm_marker_gene_dropdown,
      endoderm_featureName_series_input,
      endoderm_featureList_series_textarea,
    ]
  ),
  body=True,
  style = {'height':'100vh'},
  id = 'endoderm_series_controls'
)

endoderm_umap_controls = dbc.Card(
  [endoderm_gene_dropdown, endoderm_markerSize_dropdown],
  id = "endoderm_umap_controls",
  body=True,
  style = {'height':'86vh'},
)

endoderm_umap_tab = dbc.Tab(
  [
    dbc.Row([
      dbc.Col(
        endoderm_umap_controls,
        width = 2
      ),
      dbc.Col(
        [
          dbc.Row(
            [
              dbc.Col(
                dcc.Graph(figure=celltytpeUmapPlaceholder, id="endoderm_ctp_umap",style={'height': "80vh", 'width': '43vw'})
              ),
              dbc.Col(
                [
                  dbc.Row(
                    [
                      dbc.Col(
                        [
                          dcc.Graph(figure=stageGeneUmapPlaceholder, id="endoderm_gene_umap", style={'height': "80vh", 'width': '31vw'}),
                        ]
                      )
                    ]
                  )
                ]
              ),
            ], align = "center", style={'height': '80vh'}
          ), 
          endoderm_stage_slider
        ]
      )
    ])
  ], 
  label="UMAP",
  tab_id = "endoderm_umap_tab"
)

endoderm_series_tab = dbc.Tab(
  [
    dbc.Row(
      [
        dbc.Col(
          [
            endoderm_series_controls
          ], 
          width=2
        ),
        dbc.Col([
          dbc.Row([
            dcc.Graph(figure=celltypeUmapSeriesPlaceholder, id="endoderm_ctp_series")
          ]),
          dbc.Row([html.Img(id = 'endoderm_plotFeatureSeries_img', style = {'width': '80vw', 'margin-left':'1vw'})]),
        ], width=10)
      ]
    )
  ],
  label = "Plot Feature(Series)",
  tab_id = "endoderm_series_tab"
)

endoderm_tabs = dbc.Card(
  dbc.Tabs(
    [endoderm_umap_tab, endoderm_series_tab],
    active_tab = "endoderm_umap_tab",  
    id = "endoderm_tabs",
  )
)

layout = dbc.Container(
  [
    html.Div(id='endoderm_container'),
    endoderm_tabs
  ],
  fluid=True,
  className="dbc",
)