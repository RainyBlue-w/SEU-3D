import pandas as pd

df = pd.read_csv('/data1/share/omics-viewer/spatial/celltype_cmap.csv')
ctp_cmap = dict(zip(df['celltype'], df['color']))
ctp_cmap_reik = dict(zip(df['celltype'].str.replace(' ', '_'), df['color']))