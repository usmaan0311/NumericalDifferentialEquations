import numpy as np, pandas as pd, glob, matplotlib.pyplot as plt

files=glob.glob("*.dat")

dfs=pd.read_csv(files[0], sep='\s+', header=None,comment='#').dropna()

t,x = dfs.iloc[:,0], dfs.iloc[:,1]

plt.plot(t,x,label="Solution")
plt.legend()
plt.show()
