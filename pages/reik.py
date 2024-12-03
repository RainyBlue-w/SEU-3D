#!/usr/bin/env python
# coding: utf-8

import dash
dash.register_page(__name__)

# In[]: env
from dash import dcc, html, dash_table
from plotnine import *
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
import math
import matplotlib 
from dash_extensions.enrich import html
from callbacks.reik import *
matplotlib.use('agg')

# In[]: data

reik_dropdown_showMarker_celltype = html.Div(
  [
    dbc.Label("CellType"),
    dcc.Dropdown(
      id="reik_dropdown_showMarker_celltype",
      clearable=True,
      searchable=True,
    ),
  ],
  className="mb-4",
)

reik_datatable_showMarker_marker = html.Div(
  [
    dbc.Label("Marker Genes"),
    dash_table.DataTable(
      page_action="native",
      page_current= 0,
      page_size= 21,
      active_cell = {'row': 0, 'column': 0, 'column_id': 'Name'},
      filter_action='custom',
      filter_query='',
      filter_options={"placeholder_text": "filter..."},
      style_table={'width': '210px'},
      style_cell={'textAlign': 'center'},
      style_cell_conditional=[
        {'if': {'column_id': 'Name'},
         'width': '60%'},
        {'if': {'column_id': 'Score'},
         'width': '40%'},
      ],
      style_header={
        'fontWeight': 'bold'
      },
      id = "reik_showMarker_markerGenes"
    )
  ],
)

reik_dropdown_plotFeature_gene = html.Div(
  [
    dbc.Label("Gene to Display"),
    dcc.Dropdown(
      getStageAdata(exp_data, 'E7.5').var_names,
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
    dbc.Label("Regulon to Display"),
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
    dbc.Label("Regulon's Target Gene"),
    dash_table.DataTable(sort_action="native", page_action='native',
                         page_current= 0, page_size= 10,
                         id = "reik_plotFeature_regulonTargetGene",
                         style_cell={
                           'text-align': 'center',
                         },
                         filter_action="native", 
                         filter_options={"placeholder_text": "filter gene..."}
                        ),
  ]
)

reik_celltypeGraph = dcc.Graph(figure={}, id="reik_plotFeature_celltypeGraph")

reik_featureGraphs = dbc.Row([
  dbc.Col([
    dcc.Graph(figure={}, id='reik_plotFeature_geneGraph')
  ]),
  dbc.Col([
    dcc.Graph(figure={}, id='reik_plotFeature_regulonGraph')
  ])
])

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
          reik_table_regulonTargetGene
        ],
        body=True,
        id = 'reik_plotFeature_controlPanel'
      )
    ], width="auto"),
    # fig panel
    dbc.Col([
      dbc.Row([
        reik_celltypeGraph,
      ]),
      reik_featureGraphs,
      dbc.Row([
        dcc.Slider(id='reik_plotFeature_stageSlider',
            min = math.floor(num_stages[0]),
            max = math.ceil(num_stages[-1]),
            value=num_stages[0], step=None,
            marks= dict(zip(
              num_stages,
              stage_list))
          ),
      ], style={'margin-top': '20px'})
    ], width=9)
  ])],
  label = "Plot Feature",
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

# show marker
reik_tab_showMarker = dbc.Tab(
  [
    dbc.Row([
      dbc.Col([
        dbc.Card(
        [
          reik_dropdown_showMarker_celltype,
          reik_datatable_showMarker_marker
        ],
        style={"height": "128vh"},
        body=True,
        id = 'reik_showMarker_controlPanel'
      )
      ], width="auto"),
      dbc.Col([
        dbc.Row([
           dbc.Col([
              dcc.Graph(figure={}, id='reik_showMarker_geneGraph')
           ], width={"offset": 3})
        ]),
        dbc.Row([
           dcc.Graph(figure={}, id='reik_showMarker_markerGraph')
        ]),
        dbc.Row([
          dcc.Slider(id='reik_showMarker_stageSlider',
              min = math.floor(num_stages[0]),
              max = math.ceil(num_stages[-1]),
              value=num_stages[0], step=None,
              marks= dict(zip(
                num_stages,
                stage_list))
            ),
        ])
      ])
    ])
  ],
  label = "Show Marker",
  tab_id = "reik_tab_showMarker"
)

# In[]: app/all-layout

reik_tabs = dbc.Tabs(
  [reik_tab_plotFeature, reik_tab_showMarker, reik_tab_plotFeatureSeries],
  active_tab = "reik_tab_plotFeature",  
  id = "reik_tabs",
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