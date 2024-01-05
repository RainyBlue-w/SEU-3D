from sqlalchemy import create_engine
import os
import pandas as pd
import anndata

passwd = 'wc207279'

engine = create_engine("mysql+pymysql://wuc::207279@localhost/dash_spatial", echo=True)

connect = engine.connect()



