import pickle
import os
import scanpy as sc
import pandas as pd
import plotly.express as px

def getStageGeneExp(celltypeUmap, adata, geneList, pkl_path):
    """
        计算每个时期不同基因的表达情况
        
        Params:
        celltypeUmap(dict):键为stage，值为dataFrame，记录了celltype和对应umap坐标
        adata(Anndata):单细胞Anndata对象
        geneList(list):基因列表
        pkl_path(str):计算结果序列化数据保存位置

        Returns:
        data(dict):二重字典，第一层key为stage，第二层key为gene，值为dataframe，对应基因表达和umap坐标        
    """
    if os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        data = {}
        for stage in celltypeUmap:
            data[stage] = {}
            df = celltypeUmap[stage].loc[:, ['umapX', 'umapY']]
            stageX = adata[adata.obs['stage'] == stage].X
            for index, geneName in enumerate(geneList):
                data[stage][geneName] = df.copy()
                data[stage][geneName]['geneExp'] = list(stageX[:, index])
        dumpPkl(data, pkl_path)
        return data

def getCelltypeUmapSeriesFig(celltypeUmap, cellColor, stage):
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
    data = px.scatter(celltypeUmap[stage], x="umapX", y="umapY", color="celltype",
                      color_discrete_map=cellColor)
    data.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        yaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        legend=dict(
              title='',
              orientation='h',
              yanchor='middle',
              xanchor='right',
              y=0.5,
              x=3
        )
    )
    return data

def getGeneUmapSeriesFig(stageGeneExp, stage, gene):
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
        (0.00, "#eeeeee"),
        (0.05, "#eeeeee"),
        (1.00, "#225EA8")
    ]
    data = px.scatter(stageGeneExp[stage][gene], x="umapX", y="umapY", color="geneExp", color_continuous_scale=cmap)
    data.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=True, showline=True, showticklabels=True, title='umap1', linecolor='black', dtick=5),
        yaxis=dict(showgrid=True, showline=True, showticklabels=True, title='umap2', linecolor='black', dtick=5),
        coloraxis_colorbar_title=gene
    )
    return data

def getGeneUmapFig(stageGeneExp, stage, gene):
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
        (0.00, "#eeeeee"),
        (0.05, "#eeeeee"),
        (1.00, "#225EA8")
    ]
    data = px.scatter(stageGeneExp[stage][gene], x="umapX", y="umapY", color="geneExp", color_continuous_scale=cmap)
    data.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        yaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        coloraxis_colorbar_title=gene
    )
    return data

def getCelltypeUmapFig(celltypeUmap, cellColor, stage):
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
    data = px.scatter(celltypeUmap[stage], x="umapX", y="umapY", color="celltype",
                      color_discrete_map=cellColor)
    data.update_layout(
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        yaxis=dict(showgrid=False, showline=False, showticklabels=False, title=''),
        legend=dict(title='')
    )
    return data

def getCellTypeUmap(adata, pkl_path):
    """
        根据adata获取celltype umap坐标

        Params：
        adata(Anndata):单细胞数据
        pkl_path(str):序列化celltypeUmap坐标路径

        Returns：
        data(dict):键为stage，值为dataFrame，列celltype， umapX， umapY
    """
    if os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        data = {}
        for stage, group in adata.obs.groupby('stage'):
            data[stage] = pd.DataFrame(group[['celltype', 'umapX', 'umapY']])
        dumpPkl(data, pkl_path)
        return data

def getCellTypeColor(adata, pkl_path):
     """
        从adata中获取细胞与颜色的映射

        Params:
        adata(Anndata):单细胞数据
        pkl_path(str):序列化后celltypeColor存储路径

        Return:
        data(dict):细胞类型与颜色映射字典
     """
     if os.path.exists(pkl_path):
        return loadPkl(pkl_path)
     else:
        celltype = adata.obs['celltype']
        cellcolor = adata.obs['colour']
        data = {cell:color for cell,color in zip(celltype, cellcolor)}
        dumpPkl(data, pkl_path)
        return data
    
def getStageDict(adata, pkl_path):
    """
        根据adata获取stage

        Params:
        adata(Anndata):单细胞数据
        pkl_path(str):序列化后stageDict存储路径

        Returns:
        data(dict):胚胎发育阶段及对应数字索引
    """
    if os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        data = adata.obs['stage'].unique().tolist()
        data.sort()
        data = {i:data[i] for i in range(len(data))}
        dumpPkl(data, pkl_path)
        return data

def getGeneList(adata, pkl_path):
    """
        根据adata获取gene列表

        Params:
        adata(Anndata):单细胞数据
        pkl_path(str):序列化后geneList存储路径

        Returns:
        data(list):基因列表
    """
    if os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        data = adata.var['features'].unique().tolist()
        dumpPkl(data, pkl_path)
        return data

def loadAtlasAllData(h5ad_path, pkl_path):
    """
        加载atlasAll.h5ad数据

        Params:
        h5ad_path(str):h5ad格式数据路径
        pkl_path(str):序列化后Anndata对象存储路径

        Returns:
        adata(Anndata):单细胞数据Anndata对象
    """
    if os.path.exists(pkl_path):
        return loadPkl(pkl_path)
    else:
        adata = sc.read_h5ad(h5ad_path)
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