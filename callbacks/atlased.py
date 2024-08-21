import diskcache
from utils.atlased import *
from dash.exceptions import PreventUpdate
import re
import dash
from dash import callback, Output, Input, State, DiskcacheManager, dcc

background_callback_manager = DiskcacheManager(diskcache.Cache("/data1/share/omics-viewer/atlas-Endoderm/cache"))

# 设置根路径
root_path = '/data1/share/omics-viewer/atlas-Endoderm/'
root_pkl = '/data1/share/omics-viewer/atlas-Endoderm/pickle/'

# 设置文件名称
fileName = 'endodermCombined'
h5ad_fileName = fileName+'.h5ad'
h5ad_pklName = fileName+'_h5ad.pkl'
gene_list_pklName = fileName+'_geneList.pkl'
stage_dict_pklName = fileName+'_stageDict.pkl'
cell_color_pklName = fileName+'_cellColor.pkl'
celltype_umap_pklName = fileName+'_celltypeUmap.pkl'
# stage_geneExp_pklName = fileName+'_stageGeneExp.pkl'

# 设置文件路径
h5ad_path = root_path+h5ad_fileName
h5ad_pkl = root_pkl+h5ad_pklName
gene_list_pkl = root_pkl+gene_list_pklName
stage_dict_pkl = root_pkl+stage_dict_pklName
cell_color_pkl = root_pkl+cell_color_pklName
celltype_umap_pkl = root_pkl+celltype_umap_pklName
# stage_geneExp_pkl = root_pkl+stage_geneExp_pklName

# 加载数据
endodermData = loadData(h5ad_path, h5ad_pkl)
geneIndexDict = getGeneIndexDict(endodermData)
cellColor = getCellTypeColor(endodermData, cell_color_pkl)
geneList = getGeneList(endodermData, gene_list_pkl)
stageDict = getStageDict(endodermData, stage_dict_pkl)
stageValues = list(stageDict.values())
celltypeUmap = getCellTypeUmap(endodermData, celltype_umap_pkl)
featureName = {'celltype':'celltype', 'stage':'orig.stage', 'endo_gutCluster':'endo_gutCluster'}
# stageGeneExp = getStageGeneExp(celltypeUmap, atlasAllData, geneList, stage_geneExp_pkl)

# 初始图像占位
celltytpeUmapPlaceholder = getCelltypeUmapFig(celltypeUmap, cellColor, stageDict[2])
celltypeUmapSeriesPlaceholder = getCelltypeUmapSeriesFig(celltypeUmap, cellColor, stageDict[3])
stageGeneUmapPlaceholder = getGeneUmapFig(celltypeUmap, endodermData, stageDict[2], geneList[3], geneIndexDict)

# 回调函数

# 点击标记基因，更新feature列表
@callback(
  Output('endoderm_featureList_series_textarea', 'value', allow_duplicate=True),
  Input('endoderm_marker_gene_dropdown', 'value'),
  State('endoderm_featureList_series_textarea', 'value'),
  prevent_initial_call = True,
)
def update_featureList(gene, gene_list):
  if gene:
    if gene_list:
      if gene not in gene_list:
        gene_list += " "+gene
    else:
      gene_list = gene
  return gene_list 

# 根据套索内容获取标记基因
@callback(
  Output('endoderm_marker_gene_dropdown', 'options'),
  Output('endoderm_marker_gene_dropdown', 'value'),
  Output('endoderm_featureList_series_textarea', 'value', allow_duplicate=True),
  Input('endoderm_marker_gene_button','n_clicks'),
  State('endoderm_ctp_series', 'selectedData'),
  State('endoderm_stage_series_dropdown', 'value'),
  prevent_initial_call = True,
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('endoderm_marker_gene_dropdown', 'disabled', allow_duplicate=True),True, False),
    (Output('endoderm_marker_gene_button', 'disabled', allow_duplicate=True), True, False),
    (Output('endoderm_marker_gene_button', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('endoderm_featureList_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('endoderm_featureList_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('endoderm_featureList_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('endoderm_featureList_series_textarea', 'disabled', allow_duplicate=True),True, False),
    (Output('endoderm_featureName_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('endoderm_featureName_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('endoderm_featureName_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('endoderm_featureName_series_input', 'disabled', allow_duplicate=True),True, False),
  ],
)
def obtain_marker_gene(click, selectedData, stage):
  if click and selectedData:
     gene_list = getMarkerGenes(endodermData, celltypeUmap[stage], selectedData['points'])
     return gene_list, None, None
  else:
     raise PreventUpdate

@callback(
  Output('endoderm_plotFeatureSeries_img', 'src', allow_duplicate=True),
  Input('endoderm_featureList_series_inputButton', 'n_clicks'),
  State('endoderm_featureList_series_textarea', 'value'),
  State('endoderm_stage_series_dropdown', 'value'),
  State('endoderm_markerSize_dropdown', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('endoderm_featureList_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('endoderm_featureList_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('endoderm_featureList_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('endoderm_featureList_series_textarea', 'disabled', allow_duplicate=True),True, False),
    (Output('endoderm_featureName_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('endoderm_featureName_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('endoderm_featureName_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('endoderm_featureName_series_input', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_endoderm_plotFeature_graphSeries_pattern(click, names, stage, markerSize):
  if names is None:
    raise PreventUpdate
  if click:
    names = re.split(", |,| |\n|\'|\"|#|_|%|$|@|\(|\)|\||^|&", names)
    geneSet = set(geneList)
    genes = []
    for name in names:
       if name in geneSet:
            genes.append(name)   
    return getGeneUmapSeriesImg(celltypeUmap, endodermData, stage, genes, geneIndexDict, dot_size=2*markerSize)
  else:
    raise PreventUpdate

@callback(
  Output('endoderm_plotFeatureSeries_img', 'src', allow_duplicate=True),
  Input('endoderm_featureName_series_inputButton', 'n_clicks'),
  State('endoderm_featureName_series_input', 'value'),
  State('endoderm_stage_series_dropdown', 'value'),
  State('endoderm_markerSize_dropdown', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('endoderm_featureList_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('endoderm_featureList_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('endoderm_featureList_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('endoderm_featureList_series_textarea', 'disabled', allow_duplicate=True),True, False),
    (Output('endoderm_featureName_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('endoderm_featureName_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('endoderm_featureName_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('endoderm_featureName_series_input', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_endoderm_plotFeature_graphSeries_pattern(click, pattern, stage, markerSize):
  if pattern is None:
    raise PreventUpdate
  if click:
    pattern = pattern.lower()
    genes = []
    for gene in geneList:
       if pattern in gene.lower():
          genes.append(gene)
    return getGeneUmapSeriesImg(celltypeUmap, endodermData, stage, genes, geneIndexDict, dot_size=2*markerSize)
  else:
    raise PreventUpdate

@callback(
  Output('endoderm_geneNumber_series_text', 'children'),
  Input('endoderm_featureName_series_input', 'value'),
)
def update_endoderm_geneNumber_series(pattern):
  if not pattern:
     return '0 genes'  
  pattern = pattern.lower()
  counts = sum(1 for gene in geneList if pattern in gene.lower())
  str = f'{counts} genes'
  return str

@callback(
    Output('endoderm_ctp_series', 'figure'),
    Input('endoderm_stage_series_dropdown', 'value'),
    Input('endoderm_markerSize_dropdown', 'value'),
    Input('endoderm_featureType_series_dropdown', 'value')
)
def update_endoderm_series_ctp_umap(stage, marker_size, feature):
    return getCelltypeUmapSeriesFig(celltypeUmap, cellColor, stage, markerSize=marker_size, feature=featureName[feature])

@callback(
    Output('endoderm_gene_umap', 'figure'),
    Input('endoderm_gene_dropdown', 'value'),
    Input('endoderm_stage_slider', 'value'),
    Input('endoderm_markerSize_dropdown', 'value')
)
def update_endoderm_gene_umap(gene_name, stageIndex, marker_size):
    return getGeneUmapFig(celltypeUmap, endodermData, stageDict[stageIndex], gene_name, geneIndexDict, markerSize=marker_size)

@callback(
    Output('endoderm_ctp_umap', 'figure'),
    Input('endoderm_stage_slider', 'value'),
    Input('endoderm_markerSize_dropdown', 'value')
)
def update_endoderm_ctp_umap(stageIndex, marker_size):
    return getCelltypeUmapFig(celltypeUmap, cellColor, stageDict[stageIndex], markerSize=marker_size)