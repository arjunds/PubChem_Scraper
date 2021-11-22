import sqlite3
import pandas as pd
from ms2binner import bin

num = 10000
query = "select assays.*, spectra.* from assays as assays inner join spectra on assays.`Spectrum Index`=spectra.`index` order by random() limit " + str(num)

con = sqlite3.connect("assay_spectra.db")

df = pd.read_sql(query, con)

X,bins,scans = bin.bin_sql(df)
