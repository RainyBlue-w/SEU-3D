import pickle
import os
import scanpy as sc
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import numpy as np
from io import BytesIO
from matplotlib.colors import LinearSegmentedColormap
import base64
import anndata

def getMarkerGenes(adata, celltypeUmap_df, points):
    """
        根据套索结果计算差异基因       
    """
    pointIndex = {str(d['x'])+' '+str(d['y']) for d in points}
    allIndex = set(celltypeUmap_df.index)
    leftIndex = allIndex-pointIndex    
    points_cellid = celltypeUmap_df.loc[list(pointIndex)]['cell_id']
    left_cellid = celltypeUmap_df.loc[list(leftIndex)]['cell_id']
    selectedData = adata[points_cellid].to_memory().copy()
    leftData = adata[left_cellid].to_memory().copy()
    selectedData.obs = selectedData.obs.iloc[:, 0:0]
    selectedData.obs['group'] = 'selected'
    leftData.obs = leftData.obs.iloc[:, 0:0]
    leftData.obs['group'] = 'others'
    combined_ad = anndata.concat([selectedData, leftData], axis=0, join='outer')
    sc.tl.rank_genes_groups(combined_ad, 'group', groups=['selected'], reference='others', method='t-test', n_genes=25)
    gene_list = combined_ad.uns['rank_genes_groups']['names']['selected'].tolist()
    return gene_list

def getGeneUmapSeriesImg(celltypeUmap, adata, stage, genes, geneIndexDict, figsize=(6.4,5.2), n_cols=3, dot_size=5, umap1='stagedUmap1', umap2='stagedUmap2'):
    color_exp = [
        (0.00, "#BFBFBF"),
        (0.05, "#BFBFBF"),
        (0.75, "#225EA8"),
        (1.00, "#000000")
    ]
    cmp = LinearSegmentedColormap.from_list('custom_exp', color_exp)
    n_rows = int(np.ceil(len(genes) / n_cols))
    figsize = (figsize[0]*n_cols, figsize[1]*n_rows)
    fig = plt.figure(figsize=figsize, dpi=200)
    i = 1
    for gene in genes:
        df = getStageGeneExp(celltypeUmap, adata, stage, gene, geneIndexDict)
        ax = plt.subplot(n_rows, n_cols, i)
        i = i + 1
        umapX = df[umap1]
        umapY = df[umap2]
        geneExp = df['geneExp']
        plt.scatter(umapX, umapY, c=geneExp, cmap=cmp, s=dot_size, vmin = 0, vmax = geneExp.max())
        plt.colorbar(label='')
        plt.xlabel('UMAP-1')
        plt.ylabel('UMAP-2')
        plt.title(gene, fontsize=16)
       
    plt.tight_layout()
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    image_png = buffer.getvalue()
    graph = base64.b64encode(image_png).decode('utf-8')
    plt.close()
    return 'data:image/png;base64,{}'.format(graph)

def getStageGeneExp(celltypeUmap, adata, stage, geneName, geneIndexDict, umap1='stagedUmap1', umap2='stagedUmap2'):
    """
        计算每个时期不同基因的表达情况
        
        Params:
        celltypeUmap(dict):键为stage，值为dataFrame，记录了celltype和对应umap坐标
        adata(Anndata):单细胞Anndata对象
        stage(str):发育时期
        genename(str):基因
        geneIndexDict(dict):key->geneName, value->index

        Returns:
        data(dataframe):对应基因表达和umap坐标        
    """
    data = celltypeUmap[stage].loc[:, [umap1, umap2]]
    stageX = adata[adata.obs['stage'] == stage].X
    data['geneExp'] = stageX[:, geneIndexDict[geneName]].toarray()
    return data.sort_values(by='geneExp')

# def getStageGeneExp(celltypeUmap, adata, geneList, pkl_path):
#     """
#         计算每个时期不同基因的表达情况
        
#         Params:
#         celltypeUmap(dict):键为stage，值为dataFrame，记录了celltype和对应umap坐标
#         adata(Anndata):单细胞Anndata对象
#         geneList(list):基因列表
#         pkl_path(str):计算结果序列化数据保存位置

#         Returns:
#         data(dict):二重字典，第一层key为stage，第二层key为gene，值为dataframe，对应基因表达和umap坐标        
#     """
#     if os.path.exists(pkl_path):
#         return loadPkl(pkl_path)
#     else:
#         data = {}
#         for stage in celltypeUmap:
#             data[stage] = {}
#             df = celltypeUmap[stage].loc[:, ['umapX', 'umapY']]
#             stageX = adata[adata.obs['stage'] == stage].X
#             for index, geneName in enumerate(geneList):
#                 data[stage][geneName] = df.copy()
#                 # data[stage][geneName]['geneExp'] = list(stageX[:, index])
#                 data[stage][geneName]['geneExp'] = stageX[:, index].toarray()
#         dumpPkl(data, pkl_path)
#         return data

def getCelltypeUmapSeriesFig(celltypeUmap, cellColor, stage, umap1='stagedUmap1', umap2='stagedUmap2', markerSize=5, feature='celltype'):
    """
        绘制series tab细胞umap图

        Params:
        celltypeUmap(dict):发育阶段对应的celltype和umap坐标
        cellColor(dict):细胞对应颜色
        stage(str):绘制某个时期的细胞umap
        width(number):绘制图像的宽度

        Returns：
        data(Figure):绘制的细胞umap图
    """
    sorted_df = celltypeUmap[stage].sort_values(feature)
    sorted_df[feature] = pd.Categorical(sorted_df[feature], categories=list(sorted_df[feature].unique()))
    if feature!='celltype':
        cellColor = {'NA':'#D3D3D3'}
    data = px.scatter(sorted_df, x=umap1, y=umap2, color=feature,
                    color_discrete_map=cellColor)
    data.update_traces(marker=dict(size=markerSize))
    data.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        yaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        legend=dict(
              title='',
              itemsizing='constant',
              orientation='h',
              yanchor='middle',
              xanchor='left',
              y=0.5,
              x=1
        ),
        width=1160,
        height=580
    )
    return data

def getGeneUmapFig(celltypeUmap, adata, stage, gene, geneIndexDict, umap1='stagedUmap1', umap2='stagedUmap2', markerSize=5):
    """
        绘制基因umap图

        Params:
        stageGeneExp(dict):二重字典，第一层key为stage，第二层key为gene，值为dataframe，对应基因表达和umap坐标
        stage(str):绘制某个时期的基因umap
        gene(str):要绘制的基因名称

        Returns：
        data(Figure):绘制的基因umap图
    """
    cmap = [
        (0.00, "#BFBFBF"),
        (0.05, "#BFBFBF"),
        (0.75, "#225EA8"),
        (1.00, "#000000")
    ]
    df = getStageGeneExp(celltypeUmap, adata, stage, gene, geneIndexDict)
    data = px.scatter(df, x=umap1, y=umap2, color="geneExp", color_continuous_scale=cmap)
    data.update_traces(marker=dict(size=markerSize))
    data.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        yaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        coloraxis_colorbar_title=gene
    )
    return data

def getCelltypeUmapFig(celltypeUmap, cellColor, stage, umap1='stagedUmap1', umap2='stagedUmap2', markerSize=5):
    """
        绘制细胞umap图

        Params:
        celltypeUmap(dict):发育阶段对应的celltype和umap坐标
        cellColor(dict):细胞对应颜色
        stage(str):绘制某个时期的细胞umap
        width(number):绘制图像的宽度
        Returns：
        data(Figure):绘制的细胞umap图
    """
    sorted_df = celltypeUmap[stage].sort_values('celltype')
    sorted_df['celltype'] = pd.Categorical(sorted_df['celltype'], categories=list(sorted_df['celltype'].unique()))
    data = px.scatter(sorted_df, x=umap1, y=umap2, color="celltype",
                      color_discrete_map=cellColor)
    data.update_traces(marker=dict(size=markerSize))
    data.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        yaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        legend=dict(title='', itemsizing='constant')
    )
    return data

def getGeneIndexDict(adata):
    """
        获取gene对应索引映射

        Params:
        adata(anndata):单细胞Anndata对象

        Returns：
        data(dict):key->geneName, value->index
    """
    data = {}
    index = 0
    for geneName in adata.var.index:
        data[geneName] = index
        index+=1
    return data

def getCellTypeUmap(adata, pkl_path=None):
    """
        根据adata获取celltype umap坐标

        Params：
        adata(Anndata):单细胞数据
        pkl_path(str):序列化celltypeUmap坐标路径

        Returns：
        data(dict):键为stage，值为dataFrame，列celltype， umapX， umapY， stagedUMAP1， stagedUMAP2
    """
    if pkl_path and os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        data = {}
        for stage in adata.obs['stage'].unique():
            obs = adata.obs
            data[stage] = obs[obs['stage']==stage][['celltype', 'umapX', 'umapY', 'stagedUmap1', 'stagedUmap2', 'orig.stage', 'endo_gutCluster']]
            data[stage]['cell_id'] = data[stage].index
            data[stage].index = data[stage].apply(lambda row: f"{row['stagedUmap1']} {row['stagedUmap2']}", axis=1)
            uniqueType = data[stage]['celltype'].unique().tolist()
            data[stage]['celltype'] = pd.Categorical(data[stage]['celltype'], categories=uniqueType)
        if pkl_path:
            dumpPkl(data, pkl_path)
        return data

def getCellTypeColor(adata, pkl_path=None):
     """
        从adata中获取细胞与颜色的映射

        Params:
        adata(Anndata):单细胞数据
        pkl_path(str):序列化后celltypeColor存储路径

        Return:
        data(dict):细胞类型与颜色映射字典
     """
     if pkl_path and os.path.exists(pkl_path):
        return loadPkl(pkl_path)
     else:
        celltype = adata.obs['celltype']
        cellcolor = adata.obs['colour']
        data = {cell:color for cell,color in zip(celltype, cellcolor)}
        if pkl_path:
            dumpPkl(data, pkl_path)
        return data
    
def getStageDict(adata, pkl_path=None):
    """
        根据adata获取stage

        Params:
        adata(Anndata):单细胞数据
        pkl_path(str):序列化后stageDict存储路径

        Returns:
        data(dict):胚胎发育阶段及对应数字索引
    """
    if pkl_path and os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        data = adata.obs['stage'].unique().tolist()
        data.sort()
        data = {i:data[i] for i in range(len(data))}
        if pkl_path:
            dumpPkl(data, pkl_path)
        return data

def getGeneList(adata, pkl_path=None):
    """
        根据adata获取gene列表

        Params:
        adata(Anndata):单细胞数据
        pkl_path(str):序列化后geneList存储路径

        Returns:
        data(list):基因列表
    """
    if pkl_path and os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        data = adata.var['gene'].unique().tolist()
        if pkl_path:
            dumpPkl(data, pkl_path)
        return data

def loadData(h5ad_path, pkl_path=None):
    """
        加载atlasAll.h5ad数据

        Params:
        h5ad_path(str):h5ad格式数据路径
        pkl_path(str):序列化后Anndata对象存储路径

        Returns:
        adata(Anndata):单细胞数据Anndata对象
    """
    if not pkl_path:
        return sc.read_h5ad(h5ad_path, backed='r')
    if pkl_path and os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        adata = sc.read_h5ad(h5ad_path)
        # 根据stage计算新的UMAP坐标
        stageList = adata.obs['stage'].unique()
        stagedUmap1 = {}
        stagedUmap2 = {}
        for stage in stageList:
            subset = adata[adata.obs['stage']==stage]
            sc.tl.pca(subset, svd_solver='arpack')
            if subset.obs['sequencing.batch'].unique().size<2:
                sc.pp.neighbors(subset, n_neighbors=10, n_pcs=40)
            else:
                sc.external.pp.bbknn(subset, batch_key='sequencing.batch')
            # if '-' in stage:
            #     sc.external.pp.bbknn(subset, batch_key='orig.stage')
            sc.tl.umap(subset)
            for i, name in enumerate(subset.obs.index):
                stagedUmap1[name] = subset.obsm['X_umap'][i][0]
                stagedUmap2[name] = subset.obsm['X_umap'][i][1]
        sUmap1 = []
        sUmap2 = []
        for cellid in adata.obs.index:
            sUmap1.append(stagedUmap1[cellid])
            sUmap2.append(stagedUmap2[cellid])
        adata.obs['stagedUmap1'] = sUmap1
        adata.obs['stagedUmap2'] = sUmap2
        # adata = adata.raw.to_adata()
        if pkl_path:
            dumpPkl(adata, pkl_path)
        return adata

def loadPkl(pkl_path):
    """
        加载pkl序列化数据

        Params:
        pkl_path(str):加载序列化对象的存储路径

        Returns:
        data(Object):加载后的数据对象
    """
    with open(pkl_path, 'rb') as file:
            return pickle.load(file)
    
def dumpPkl(data, pkl_path):
     """
        写入pkl序列化数据
        
        Params:
        data(Object):需要序列化的对象
        pkl_path(str):序列化对象的存储路径

        Returns:
        None
     """
     with open(pkl_path, 'wb') as file:
            pickle.dump(data, file)