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
from callbacks.atlasall import *

background_callback_manager = DiskcacheManager(diskcache.Cache('/data1/share/omics-viewer/atlas-ALL/cache'))

atlasall_gene_dropdown = html.Div(
  [
    dbc.Label("Gene to Display"),
    dcc.Dropdown(
      geneList,
      geneList[3],
      id="atlasall_gene_dropdown",
      clearable=False,
      searchable=True,
    ),
  ], className="mb-4",
)

atlasall_stage_slider = html.Div(
  dcc.Slider(
    0, len(stageDict)-1,
    step = None,
    marks=stageDict,
    value=2,
    id="atlasall_stage_slider"
  ), style={'width': '97%'}
)

atlasall_featureType_series_dropdown = html.Div(
  [
    dbc.Label("Feature Type"),
    dcc.Dropdown(
      ['Gene'],
      'Gene',
      id="atlasall_featureType_series_dropdown",
      clearable=False,
      searchable=True,
    ),
  ], className="mb-4",
)

atlasall_stage_series_dropdown = html.Div(
  [
    dbc.Label("Stage"),
    dcc.Dropdown(
      stageValues,
      stageValues[3],
      id="atlasall_stage_series_dropdown",
      clearable=False,
      searchable=True,
    ),
  ], className="mb-4",
)

atlasall_featureName_series_input = html.Div(
  [
    dbc.Label("Contains"),
    dbc.InputGroup(
      [
        dmc.Grid(
          children=[
            dmc.Col(dbc.Input(id="atlasall_featureName_series_input"), span=9),
            dmc.Col(dbc.Button('Plot', id='atlasall_featureName_series_inputButton', n_clicks=0, color='primary'), span=3),
            dmc.Col(dmc.Text(id='atlasall_geneNumber_series_text', color='gray'), span=12),
          ], gutter=3
        )
      ]
    )
  ], className="mb-4",
)

atlasall_featureList_series_textarea = html.Div(
  [
    html.Div(
      dmc.Grid(
        [
          dmc.Col(dbc.Label("Name List"), span=9),
          dmc.Col(dbc.Button('Plot', id='atlasall_featureList_series_inputButton', n_clicks=0, color='primary'), span=3)
        ], gutter=3
      ), className="mb-2"
    ),
    html.Div(
      dbc.Textarea(id = "atlasall_featureList_series_textarea",
        placeholder="paste feature names here(seperated by any charactor)\ne.g.  A  B,C,  D\nE##G_H,#I@KK%%G%(U),(V)|(W)\"X\",\"Y\"^Q*I",
        rows=8, className="mb-3",)
    )
  ], className="mb-4",
)

atlasall_series_controls = dbc.Card(
  dbc.Col(
    [
      atlasall_featureType_series_dropdown,
      atlasall_stage_series_dropdown,
      atlasall_featureName_series_input,
      atlasall_featureList_series_textarea,
    ]
  ),
  body=True,
  style = {'height':'86vh'},
  id = 'atlasall_series_controls'
)

atlasall_umap_controls = dbc.Card(
  atlasall_gene_dropdown,
  id = "atlasall_umap_controls",
  body=True,
  style = {'height':'86vh'},
)

atlasall_umap_tab = dbc.Tab(
  [
    dbc.Row([
      dbc.Col(
        atlasall_umap_controls,
        width = 2
      ),
      dbc.Col(
        [
          dbc.Row(
            [
              dbc.Col(
                dcc.Graph(figure=celltytpeUmapPlaceholder, id="atlasall_ctp_umap",style={'height': "80vh", 'width': '43vw'})
              ),
              dbc.Col(
                [
                  dbc.Row(
                    [
                      dbc.Col(
                        [
                          dcc.Graph(figure=stageGeneUmapPlaceholder, id="atlasall_gene_umap", style={'height': "80vh", 'width': '31vw'}),
                        ]
                      )
                    ]
                  )
                ]
              ),
            ], align = "center", style={'height': '80vh'}
          ), 
          atlasall_stage_slider
        ]
      )
    ])
  ], 
  label="UMAP",
  tab_id = "atlasall_umap_tab"
)

atlasall_series_tab = dbc.Tab(
  [
    dbc.Row(
      [
        dbc.Col(
          [
            atlasall_series_controls
          ], 
          width=2
        ),
        dbc.Col([
          dbc.Row([
            dcc.Graph(figure=celltypeUmapSeriesPlaceholder, id="atlasall_ctp_series")
          ]),
          dbc.Row([html.Img(id = 'atlasall_plotFeatureSeries_img', style = {'width': '80vw', 'margin-left':'1vw'})]),
        ], width=10)
      ]
    )
  ],
  label = "Plot Feature(Series)",
  tab_id = "atlasall_series_tab"
)

atlasall_tabs = dbc.Card(
  dbc.Tabs(
    [atlasall_umap_tab, atlasall_series_tab],
    active_tab = "atlasall_umap_tab",  
    id = "atlasall_tabs",
  )
)

layout = dbc.Container(
  [
    html.Div(id='atlasall_container'),
    atlasall_tabs
  ],
  fluid=True,
  className="dbc",
)