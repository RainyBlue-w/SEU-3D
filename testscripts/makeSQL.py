
# %%
import os
import pandas as pd
import scanpy as sc
import pymysql
from sqlalchemy import create_engine
from sqlalchemy.types import NVARCHAR, Float, Integer

# %%
exp_data = {
  'E7.5': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.5_HC0.5_min400.h5ad"),
  'E7.75': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E7.75_HC0.5_min400.h5ad"),
  'E8.0': sc.read_h5ad("/rad/wuc/dash_data/spatial/matrix_data/E8.0_HC0.5_min400.h5ad")
}

# for stage, data in exp_data.items():

coord_data = {
  'E7.5': pd.read_csv('/rad/wuc/dash_data/spatial/spatial_coordinate.embryo_E7.5.csv', sep=' '),
  'E7.75': pd.read_csv('/rad/wuc/dash_data/spatial/spatial_coordinate.embryo_E7.75.csv', sep=' '),
  'E8.0': pd.read_csv('/rad/wuc/dash_data/spatial/spatial_coordinate.embryo_E8.0.csv', sep=' ')
}

for stage in exp_data.keys():
  exp_data[stage].obs[['x','y','z','x_flatten', 'y_flatten']] = coord_data[stage].loc[exp_data[stage].obs_names,['x','y','z','x_flatten', 'y_flatten']]

for stage in exp_data.keys():
  exp_data[stage].raw = exp_data[stage].copy()
  sc.pp.normalize_total(exp_data[stage], target_sum=1e4)
  sc.pp.log1p(exp_data[stage])
  exp_data[stage].obsm = {'X_spatial': coord_data[stage].loc[exp_data[stage].obs_names,['x','y','z']]}

for stage in exp_data.keys():
  exp_data[stage]= exp_data[stage][exp_data[stage].obs_names.sort_values()]

for k,v in exp_data.items():
  v.obs.germ_layer = [i.replace('exe-ecto', 'ectoderm') for i in v.obs.germ_layer.values]

# %% sqlalchemy
engine = create_engine('mysql+mysqlconnector://wuc:207279@localhost/dash_spatial')
# %%
for stage,data in exp_data.items():
    df = data.to_df()
    df.to_sql(
        name=f'{stage}_EXP', con=engine, index=False, 
        if_exists='replace', chunksize=1000)

# %%
