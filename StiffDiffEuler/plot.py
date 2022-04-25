import numpy as np, pandas as pd, glob, matplotlib.pyplot as plt
file=glob.glob("*.dat")

df=pd.read_csv(file[0],sep='\s+',header=None).dropna()

x,y=df.iloc[:,0], df.iloc[:,1]

plt.plot(x,y,label="solution")
plt.show()
