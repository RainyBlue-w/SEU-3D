import diskcache
from utils.atlasall import *
from dash.exceptions import PreventUpdate
import re
from dash import callback, Output, Input, State, DiskcacheManager, dcc

background_callback_manager = DiskcacheManager(diskcache.Cache("/rad/share/omics-viewer/atlas-CP/cache"))

# 设置根路径
root_path = '/rad/share/omics-viewer/atlas-ALL/'
root_pkl = '/rad/share/omics-viewer/atlas-ALL/pickle/'

# 设置文件名称
fileName = 'atlasAll'
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
atlasAllData = loadAtlasAllData(h5ad_path, h5ad_pkl)
geneIndexDict = getGeneIndexDict(atlasAllData)
cellColor = getCellTypeColor(atlasAllData, cell_color_pkl)
geneList = getGeneList(atlasAllData, gene_list_pkl)
stageDict = getStageDict(atlasAllData, stage_dict_pkl)
stageValues = list(stageDict.values())
celltypeUmap = getCellTypeUmap(atlasAllData, celltype_umap_pkl)
# stageGeneExp = getStageGeneExp(celltypeUmap, atlasAllData, geneList, stage_geneExp_pkl)

# 初始图像占位
celltytpeUmapPlaceholder = getCelltypeUmapFig(celltypeUmap, cellColor, stageDict[2])
celltypeUmapSeriesPlaceholder = getCelltypeUmapSeriesFig(celltypeUmap, cellColor, stageDict[3])
stageGeneUmapPlaceholder = getGeneUmapFig(celltypeUmap, atlasAllData, stageDict[2], geneList[3], geneIndexDict)

# 回调函数
@callback(
  Output('atlasall_plotFeatureSeries_img', 'src', allow_duplicate=True),
  Input('atlasall_featureList_series_inputButton', 'n_clicks'),
  State('atlasall_featureList_series_textarea', 'value'),
  State('atlasall_stage_series_dropdown', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('atlasall_featureList_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('atlasall_featureList_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('atlasall_featureList_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlasall_featureList_series_textarea', 'disabled', allow_duplicate=True),True, False),
    (Output('atlasall_featureName_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('atlasall_featureName_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('atlasall_featureName_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlasall_featureName_series_input', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_atlasall_plotFeature_graphSeries_pattern(click, names, stage):
  if names is None:
    raise PreventUpdate
  if click:
    names = re.split(", |,| |\n|\'|\"|#|_|%|$|@|\(|\)|\||^|&", names)
    geneSet = set(geneList)
    genes = []
    for name in names:
       if name in geneSet:
            genes.append(name)   
    return getGeneUmapSeriesImg(celltypeUmap, atlasAllData, stage, genes, geneIndexDict)
  else:
    raise PreventUpdate

@callback(
  Output('atlasall_plotFeatureSeries_img', 'src', allow_duplicate=True),
  Input('atlasall_featureName_series_inputButton', 'n_clicks'),
  State('atlasall_featureName_series_input', 'value'),
  State('atlasall_stage_series_dropdown', 'value'),
  background=True,
  manager=background_callback_manager,
  running=[
    (Output('atlasall_featureList_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('atlasall_featureList_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('atlasall_featureList_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlasall_featureList_series_textarea', 'disabled', allow_duplicate=True),True, False),
    (Output('atlasall_featureName_series_inputButton', 'children', allow_duplicate=True), 'Plot', 'Plot'),
    (Output('atlasall_featureName_series_inputButton', 'disabled', allow_duplicate=True), True, False),
    (Output('atlasall_featureName_series_inputButton', 'color', allow_duplicate=True), 'danger', 'primary'),
    (Output('atlasall_featureName_series_input', 'disabled', allow_duplicate=True),True, False),
  ],
  prevent_initial_call = True,
)
def update_atlasall_plotFeature_graphSeries_pattern(click, pattern, stage):
  if pattern is None:
    raise PreventUpdate
  if click:
    pattern = pattern.lower()
    genes = []
    for gene in geneList:
       if pattern in gene.lower():
          genes.append(gene)
    return getGeneUmapSeriesImg(celltypeUmap, atlasAllData, stage, genes, geneIndexDict)
  else:
    raise PreventUpdate

@callback(
  Output('atlasall_geneNumber_series_text', 'children'),
  Input('atlasall_featureName_series_input', 'value'),
)
def update_atlasall_geneNumber_series(pattern):
  if not pattern:
     return '0 genes'  
  pattern = pattern.lower()
  counts = sum(1 for gene in geneList if pattern in gene.lower())
  str = f'{counts} genes'
  return str

@callback(
    Output('atlasall_ctp_series', 'figure'),
    Input('atlasall_stage_series_dropdown', 'value'),
)
def update_atlasall_series_gene_umap(stage):
    return getCelltypeUmapSeriesFig(celltypeUmap, cellColor, stage)

@callback(
    Output('atlasall_gene_umap', 'figure'),
    Input('atlasall_gene_dropdown', 'value'),
    Input('atlasall_stage_slider', 'value')
)
def update_atlasall_gene_umap(gene_name, stageIndex):
    return getGeneUmapFig(celltypeUmap, atlasAllData, stageDict[stageIndex], gene_name, geneIndexDict)

@callback(
    Output('atlasall_ctp_umap', 'figure'),
    Input('atlasall_stage_slider', 'value')
)
def update_atlasall_ctp_umap(stageIndex):
    return getCelltypeUmapFig(celltypeUmap, cellColor, stageDict[stageIndex])