import os
import scanpy as sc
import pickle
import pandas as pd
import loompy as lp
import numpy as np
import re
import plotly.express as px
from plotnine import *
from io import BytesIO
import matplotlib.pyplot as plt
import matplotlib
import base64
import plotly.graph_objects as go
from sortedcontainers import SortedSet

from matplotlib.colors import LinearSegmentedColormap
color_exp = [(0.00, "#eeeeee"),
              (0.05, "#eeeeee"),
              (1.00, "#225EA8")
            ]
cmap_exp = LinearSegmentedColormap.from_list('custom_exp', color_exp)

color_auc = [(0.00, "#eeeeee"),
              (0.05, "#eeeeee"),
              (1.00, "#d53e4f")
            ]
cmap_auc = LinearSegmentedColormap.from_list('custom_auc', color_auc)

def formatStage(stage):
  """
    根据滑动条数值格式化stage
  """
  if np.remainder(stage, 0.5) == 0:
    stage = 'E%.1f' % stage
  elif np.remainder(stage, 0.25) == 0:
    stage = 'E%.2f' % stage
  return stage

def numeric_str(s):
    """
      判断一个数字是否可以转换成数值
    """
    try:
        num = float(s)
        return num
    except ValueError:
        return -100

def set_num_stages(stage_list):
  """
    预处理num_stages
  """
  num_stages = [ float(re.sub('E','',i)) for i in stage_list]
  num_stages = [int(i) if (i % 1 == 0) else i for i in num_stages ]
  return num_stages

def set_color(color_path):
  """
    设置细胞颜色
  """
  colors = pd.read_csv(color_path)
  colors['celltype'] = [re.sub(' |/', '_', i) for i in colors['celltype']]
  celltype_colors = dict(zip(colors['celltype'], colors['colour']))
  return celltype_colors

def get_file_path(file_folder_name):
  """
    根据文件或文件夹获取完整地址

    params：
    file_folder_name (str): 文件或文件夹名称

    return full_path(str) : 文件或文件夹完整路径
  """
  root = '/rad/share/omics-viewer/reik/'
  return root+file_folder_name

def load_exp_data(path, pickle_root, stage_list):
  """
    对h5ad数据做处理，并将结果序列化存储，下次可直接读取

    params：
    path(str): h5ad数据文件地址
    pickle_root（str）： 序列化数据存放路径
    stage_list(list): 发育阶段

    return exp_data(dict): key->stage; value->anadata;
  """
  picklePath = pickle_root+os.path.basename(path)+".pickle"
  if(os.path.exists(picklePath)):
    with open(picklePath, "rb") as f:
      return pickle.load(f)
  else:
    exp_data = sc.read_h5ad(path)
    exp_data = exp_data[(exp_data.obs['sample'] != 'E8.5_CRISPR_T_WT')&(exp_data.obs['sample'] != 'E8.5_CRISPR_T_KO')]
    sc.pp.normalize_total(exp_data, target_sum=1e4)
    sc.pp.log1p(exp_data)
    tmp = {}
    for stage in stage_list:
      tmp[stage] = exp_data[exp_data.obs['stage'] == stage]
    exp_data = tmp.copy()
    del(tmp)
    with open(picklePath, "wb") as f:
      pickle.dump(exp_data, f)
    return exp_data
  
def load_umap_data(umap_path, stage_list):
  """
    获取umap
  """
  umap = pd.read_csv(umap_path, index_col=0)
  tmp = {}
  for stage in stage_list:
    tmp[stage] = umap[umap.index.str.startswith(stage)]
  umap = tmp
  del(tmp)
  return umap

def load_regulon_geneset_and_auc_data(pickle_root, scenic_root, stage_list, exp_data):
  """
    加载/序列化regulon_geneset和auc_data
    params:
    pickle_root(str): 序列化后数据存储路径
    scenic_root(str): scenic路径
    return regulon_geneset, auc_data
  """
  regulon_geneset_path = pickle_root+"regulon_geneset.pickle"
  auc_data_path = pickle_root+"auc_data.pickle"
  auc_mtx = {}
  auc_data = {}
  regulon_geneset = {}
  if(os.path.exists(regulon_geneset_path)):
    with open(regulon_geneset_path, "rb") as f:
       regulon_geneset = pickle.load(f)
  else:
    for stage in stage_list:
      fpath = scenic_root+f"Reik_{stage}_pyscenic_output.loom"
      lf = lp.connect( fpath, mode='r+', validate=False )
      auc_mtx[stage] = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
      df = pd.DataFrame(lf.ra.Regulons, index = lf.ra.Gene).astype("bool")
      regulon_geneset[stage] = {}
      for regulon in df.columns:
        regulon_geneset[stage][regulon] = lf.ra.Gene[df.loc[:,regulon]].tolist()
      lf.close()
    with open(regulon_geneset_path, "wb") as f:
      pickle.dump(regulon_geneset, f)
  if(os.path.exists(auc_data_path)):
    with open(auc_data_path, "rb") as f:
       auc_data = pickle.load(f) 
  else:
    for stage in stage_list:
      auc_data[stage] = sc.AnnData(X=auc_mtx[stage], obs=exp_data[stage].obs)
    del(auc_mtx)
    with open(auc_data_path, "wb") as f:
      pickle.dump(auc_data, f)
  return regulon_geneset, auc_data

def show_feature_umap(adata, feature, embedding, cmap = None, sort=False, ascending=True, figwidth=510, figheight=300, **kws):
    """
      绘制特征umap
    """
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
    		color_continuous_scale = cmap, width=figwidth, height=figheight,
        **kws
    	)
    plot.update_yaxes(visible=False)
    plot.update_xaxes(visible=False)
    plot.update_traces(marker_size=4,
                      marker_opacity=1)
    plot.update_layout(
      margin=dict(l=0, t=0, b=0),
      plot_bgcolor = '#ffffff', 
      uirevision='constant',
      legend_itemsizing = 'constant',
      legend=dict(
            title='',
            orientation='v',
            yanchor='middle',
            xanchor='right',  # 设置图例在右侧
            y=0.5,
            x=1,  # 调整图例在横向的位置
    ),
    )
    plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],font_size = 20)) 
    return plot

def show_celltype_umap(adata, embedding, cmap, sort=False, figwidth=1150, figheight=300, **kws):
  """
    绘制细胞umap
  """
  embedding.columns = ['umap_1', 'umap_2']
  pdf = pd.concat([embedding, adata.obs[['celltype']]], axis=1)
  if sort:
    pdf = pdf.sort_values(by='celltype')
  plot = px.scatter(
  	data_frame = pdf,
    x = 'umap_1', y = 'umap_2', color = 'celltype',
    color_discrete_map = cmap, width=figwidth, height=figheight,
    **kws
  )
  plot.update_yaxes(visible=False)
  plot.update_xaxes(visible=False)
  plot.update_traces(marker_size=3, marker_opacity=1)
  plot.update_layout(
    margin=dict(l=0, t=0, b=0),
    plot_bgcolor = '#ffffff',
    title = '',
    legend=dict(
            title='CellType',
            orientation='h',
            yanchor='middle',
            xanchor='right',  # 设置图例在右侧
            y=0.5,
            x=2.7,  # 调整图例在横向的位置
    ),
    uirevision='constant',
    legend_itemsizing = 'constant'
  )
  plot.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1],font_size = 20)) 
  return plot

def show_celltype_umap_series(adata, embedding, cmap, sort=False, **kws):
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

def calculateMarker(h5ad_path, pkl_path, stage_list):
    """
    计算marker genes，之前已经计算好了就不会再计算，直接加载，否则重新计算

    params：
    h5ad_path (str): 读取的anadata文件路径。
    pkl_path (str): 序列化wilcoxon检验后结果路径。
    stage_list （list）: 胚胎发育阶段列表

    return：
    marker (dict): 键为发育阶段，值为对应时期wilcoxon检验后的结果
    """
    marker = {}
    for stage in stage_list:
        file = os.path.join(pkl_path, stage)
        if os.path.exists(file):
            with open(file, 'rb') as f:
                marker[stage] = pickle.load(f)
        else:
            adata = sc.read_h5ad(h5ad_path)
            adata = adata[(adata.obs['sample'] != 'E8.5_CRISPR_T_WT')&(adata.obs['sample'] != 'E8.5_CRISPR_T_KO')]
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            tmp = {}
            for st in stage_list:
                tmp[st] = adata[adata.obs['stage'] == st]
                filtered_data = tmp[st].obs['celltype'].value_counts()>=5
                tmp[st] = tmp[st][filtered_data[tmp[st].obs['celltype']]]
                sc.tl.rank_genes_groups(tmp[st], "celltype", method="wilcoxon")
                marker[st] = tmp[st]
                file = os.path.join(pkl_path, st)
                with open(file, 'wb') as f:
                    pickle.dump(marker[st], f)
    return marker

def get_marker_tableData(marker, stage_list, pkl_path):
    """
        获取不同阶段对应细胞类型的marker基因tabledata格式

        params:
        marker (dict): 键为发育阶段，值为对应时期wilcoxon检验后的结果
        stage_list (list): 发育时期
        pkl_path (str): 序列化后结果存储路径:/rad/zhouyb/pyProject/multi_page/reik/marker

        return df (dict): 两重字典，key1->stage, key2->celltype, value->DataFrame("marker name", "score")
    """
    file = os.path.join(pkl_path, "marker_tableData_dataFrame")
    df = {}
    if (os.path.exists(file)):
        with open(file, "rb") as f:
            df = pickle.load(f)
    else:
        for stage in stage_list:
            df[stage] = {}
            for celltype in marker[stage].obs.celltype.unique():
                data = {
                    'Name': marker[stage].uns['rank_genes_groups']['names'][celltype][:100],
                    'Score': marker[stage].uns['rank_genes_groups']['scores'][celltype][:100]
                } 
                tmp = pd.DataFrame(data)
                tmp['Score'] = tmp['Score'].round(2).astype(str).astype(float)
                df[stage][celltype] = tmp
        with open(file, "wb") as f:
            pickle.dump(df, f)
    return df

def sort_cell_by_score(marker, stage, gene):
    """
        根据指定marker gene得分对细胞进行排序
        
        params:
        marker (dict): 键为发育阶段，值为对应时期wilcoxon检验后的结果
        stage (str): 胚胎发育阶段
        gene (str): marker gene名称

        return ranked_cell_list (SortedSet): 按照marker gene得分从大到小排序后的cell结果 
    """
    cell_list = marker[stage].obs.celltype.unique()
    ranked_cell_list = SortedSet(key=lambda elem : -elem[1])
    for cell in cell_list:
        gene_list = marker[stage].uns['rank_genes_groups']['names'][cell]
        score_list = marker[stage].uns['rank_genes_groups']['scores'][cell]
        gene_index = np.where(gene_list == gene)
        score = score_list[gene_index] 
        ranked_cell_list.add((cell, score[0]))
    return ranked_cell_list

def show_features_series_matplotlib(adata, embedding, features=None, pattern=None, sort=True, ascending=True, 
                                          figsize=(6.4,4.8), n_cols=1, dot_size=4, cmap=cmap_exp, **kws):
  
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

def darken_color(hex_color, factor=0.8):
    """
      将细胞显示颜色做暗化处理，用于做小提琴图的边框颜色

      params：
      hex_color (str): 细胞颜色的十六进制值
      factor (float): 暗化程度

      return dark_hex_color (str): 暗化后的颜色值
    """
    r = int(hex_color[1:3], 16)
    g = int(hex_color[3:5], 16)
    b = int(hex_color[5:7], 16)

    r = max(0, int(r * factor))
    g = max(0, int(g * factor))
    b = max(0, int(b * factor))

    dark_hex_color = '#{:02x}{:02x}{:02x}'.format(r, g, b)

    return dark_hex_color

def show_marker_violin(marker, stage, gene, celltype_colors, w=1200, h=600):
    """
        绘制marker gene小提琴图

        params:
        marker (dict): 键为发育阶段，值为对应时期wilcoxon检验后的结果
        stage (str): 胚胎发育阶段
        gene (str): marker gene名称

        return fig (graph_objects): marker gene小提琴图 
        
    """
    ranked_cell_list = sort_cell_by_score(marker, stage, gene)
    fig = go.Figure()
    for cell in ranked_cell_list:
        cell_type = cell[0]
        df = marker[stage][marker[stage].obs.celltype==cell_type, gene].to_df()
        fig.add_trace(go.Violin(
            y=df[gene],
            name=cell_type,
            fillcolor=celltype_colors[cell_type],
            line_color=darken_color(celltype_colors[cell_type]),
            line_width=1,
            box_visible=True,
            meanline_visible=True))
        
    fig.update_layout(
        width=w,
        height=h,
        xaxis_tickangle=45,
        title={
            'text': gene,
            'x': 0.5, 
            'y': 0.9
        },
        legend={
            'x': 1,
            'y': 1,
            'orientation': 'h', 
        }
    )
    return fig

def show_features_reik_regularExp(adata, stage, odir, featureType, embedding, pattern=None, features=None, sort=False, ascending=True, dpi=100, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  if not features:
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()
  else:
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)
  pdf = pd.DataFrame(np.array(embedding), 
                      index=adata.obs_names, 
                      columns=['x', 'y'])
  features_df = adata[:, features].to_df()
  pdf = pd.concat([pdf, features_df], axis=1)
  pdf = pd.melt(pdf, id_vars = ['x', 'y'], var_name='feature', value_name='value')
  if sort:
    pdf = pdf.sort_values(by='value', ascending=ascending)
  pdf['feature'] = pdf['feature'].astype('category').values.reorder_categories(features)

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

def show_featuresCtpcounts_reik_regularExp(adata, stage, odir, featureType, embedding, cmap, pattern=None, features=None, sort=False, ascending=True, dpi=150, **kws):
  embedding = embedding.loc[adata.obs_names,:]
  if not features:
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, pattern, featureType, dpi)
    if os.path.exists( img_dir ):
      return img_dir
    features = [i  for i in adata.var_names if re.match(pattern, i)]
    features.sort()
  else:
    img_dir = '/rad/wuc/dash_data/reik/tmp/%s/%s/plotFeatureSeries_%s_%s_%s_dpi%d.png' % (odir, stage, stage, "tmp", featureType, dpi)
  ctp_counts = {}
  # for gene in ordered_features:
  for gene in features:
    df = adata[:,gene].to_df()
    df = df[df[gene] > 0]
    counts = pd.DataFrame(adata.obs['celltype'].loc[df.index].value_counts())
    counts['gene'] = gene
    counts['count'] = counts['count']/sum(counts['count'])
    ctp_counts[gene] = counts
  ctp_counts = pd.concat(ctp_counts, axis=0)
  ctp_counts['celltype'] = np.array(ctp_counts.index.to_list())[:,1]
  ctp_counts['text_y'] = ctp_counts['count'].max()/2
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