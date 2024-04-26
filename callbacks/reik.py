from utils.reik import *
import pandas as pd
from dash import Input, Output, callback, no_update, State, DiskcacheManager
from dash.exceptions import PreventUpdate
from plotnine import *
import dash_mantine_components as dmc
from dash_iconify import DashIconify
import diskcache
import dash
import re
background_callback_manager = DiskcacheManager(diskcache.Cache("/rad/share/omics-viewer/reik/cache"))

# 数据路径
h5ad_path = get_file_path('anndata.h5ad')
pickle_root = get_file_path('pickle/')
marker_pkl_path = get_file_path('marker/')
umap_path = get_file_path('all_combine_umap.csv')
scenic_root = get_file_path('scenic/')
color_path = get_file_path('celltype_color.csv')

# 发育阶段
stage_list = ['E7.5','E7.75','E8.0','E8.5','E8.75']

# 加载数据
exp_data = load_exp_data(h5ad_path, pickle_root, stage_list)
umap = load_umap_data(umap_path, stage_list)
regulon_geneset, auc_data = load_regulon_geneset_and_auc_data(pickle_root, scenic_root, stage_list, exp_data)
celltype_colors = set_color(color_path)
num_stages = set_num_stages(stage_list)
marker = calculateMarker(h5ad_path, marker_pkl_path, stage_list)
df = get_marker_tableData(marker, stage_list, marker_pkl_path)

# plot feature
@callback(
  Output('reik_dropdown_plotFeature_gene', 'options'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_gene', 'search_value')
)
def update_dropdown_options_gene(stage, search):
  if not stage:
    raise PreventUpdate
  stage = formatStage(stage)
  vars = exp_data[stage].var_names
  return vars[vars.str.lower().str.startswith(search.lower() if search else "")].sort_values()

@callback(
  Output('reik_dropdown_plotFeature_regulon', 'options'),
  Output('reik_dropdown_plotFeature_regulon', 'value'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_regulon', 'search_value'),
  State('reik_dropdown_plotFeature_regulon', 'value')
)
def update_dropdown_options_regulon(stage, search, regulon):
  if not stage:
    raise PreventUpdate
  ctx = dash.callback_context
  trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
  stage = formatStage(stage)
  vars = auc_data[stage].var_names
  value = dash.no_update if regulon in vars else vars[0]
  return vars[vars.str.lower().str.startswith(search.lower() if search else "")].sort_values(), value

@callback(
  Output('reik_plotFeature_celltypeGraph','figure'),
  Input('reik_plotFeature_stageSlider', 'value'),
)
def update_celltypeGraph(stage):
  if not stage:
    raise PreventUpdate
  stage = formatStage(stage)
  return show_celltype_umap(adata=exp_data[stage], embedding=umap[stage], cmap=celltype_colors, sort=True)

@callback(
  Output('reik_plotFeature_geneGraph','figure'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_gene', 'value')
)
def update_geneGraph(stage, gene):
  if not stage or not gene:
    raise PreventUpdate
  stage = formatStage(stage)
  return show_feature_umap(adata=exp_data[stage], feature=gene, embedding=umap[stage], sort=True)

@callback(
  Output('reik_plotFeature_regulonGraph','figure'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_regulon', 'value')
)
def update_regulonGraph(stage, regulon):
  if not stage or not regulon:
    raise PreventUpdate
  stage = formatStage(stage)
  return show_feature_umap(adata=auc_data[stage], feature=regulon, embedding=umap[stage], sort=True)

@callback(
  Output('reik_plotFeature_regulonTargetGene', 'data'),
  Output('reik_plotFeature_regulonTargetGene', 'columns'),
  Input('reik_plotFeature_stageSlider', 'value'),
  Input('reik_dropdown_plotFeature_regulon', 'value')
)
def update_regulonTargetGenes(stage, regulon):
  if not stage or not regulon:
    raise PreventUpdate
  stage = formatStage(stage)
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
  stage = formatStage(stage)
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
      adata = auc_data[stage]

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
      adata = auc_data[stage]

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

# show marker
@callback(
    Output('reik_showMarker_markerGenes', 'active_cell'),
    Input('reik_showMarker_markerGenes', "filter_query"),
)
def default_activate_cell(filter):
  """
    返回默认激活marker gene单元格
  """
  return {'row': 0, 'column': 0, 'column_id': 'Name'}

@callback(
    Output('reik_showMarker_markerGraph', 'figure'),
    Input('reik_showMarker_markerGenes', 'active_cell'),
    Input('reik_showMarker_markerGenes', 'data'),
    State('reik_showMarker_stageSlider', 'value')
)
def update_marker_markerGraph(activate_cell, data, stage):
   """
    更新marker Graph
   """
   if not data or not activate_cell:
      raise PreventUpdate
   stage = formatStage(stage)
   name = data[int(activate_cell['row'])]['Name']
   return show_marker_violin(marker, stage, name, celltype_colors, 1225, 560)


@callback(
    Output('reik_showMarker_geneGraph', 'figure'),
    Input('reik_showMarker_markerGenes', 'active_cell'),
    Input('reik_showMarker_markerGenes', 'data'),
    State('reik_showMarker_stageSlider', 'value')
)
def update_marker_geneGraph(activate_cell, data, stage):
   """
    更新marker gene Graph
   """
   if not data or not activate_cell:
      raise PreventUpdate
   stage = formatStage(stage)
   name = data[int(activate_cell['row'])]['Name']
   return show_feature_umap(adata=marker[stage], feature=name, embedding=umap[stage], sort=True, figwidth=510, figheight=300) #1020, 600

@callback(
    Output('reik_showMarker_markerGenes', 'data'),
    Output('reik_showMarker_markerGenes', 'columns'),
    Input('reik_dropdown_showMarker_celltype', 'value'),
    Input('reik_showMarker_stageSlider', 'value'),
    Input('reik_showMarker_markerGenes', "filter_query")
)
def update_marker_genes(celltype, stage, filter):
  """
    更新marker genes
  """
  if not celltype or not stage:
    raise PreventUpdate
  stage = formatStage(stage)
  if filter:
    filter = filter.replace("icontains", "scontains")
    key1 = ''
    key2 = -100
    if "&&" in filter:
      keywords = filter.split(" scontains ")
      key1 = keywords[1].split(" && ")[0]
      key2 = numeric_str(keywords[2])
    elif filter.startswith("{Name}"):
      key1 = filter.split("{Name} scontains ")[1]
    else:
      key2 = numeric_str(filter.split("{Score} scontains ")[1])
    filtered_data = df[stage][celltype][(df[stage][celltype]['Name'].str.contains('(?i)^'+key1)) & (df[stage][celltype]['Score'] > key2)]
    return filtered_data.to_dict('records'), [{"name": i, "id": i} for i in filtered_data.columns]
  return df[stage][celltype].to_dict('records'), [{"name": i, "id": i} for i in df[stage][celltype].columns]

@callback(
    Output('reik_dropdown_showMarker_celltype', 'options'),
    Output('reik_dropdown_showMarker_celltype', 'value'),
    Input('reik_showMarker_stageSlider', 'value'),
)
def update_showMarker_dropdown_celltype(stage):
  """
    更新细胞类型
  """
  if not stage:
    raise PreventUpdate
  stage = formatStage(stage)
  celltype = marker[stage].obs['celltype'].unique().sort_values()
  cell = celltype[0]
  return celltype, cell