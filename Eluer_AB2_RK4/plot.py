import pandas as pd, glob, matplotlib.pyplot as plt

files=glob.glob("*.dat")

dfs = pd.read_csv(files[1], sep='\s+', comment='#',header=None).dropna()
dfe = pd.read_csv(files[0], sep='\s+', comment='#',header=None).dropna()

t,y = dfs.iloc[:,0], dfs.iloc[:,1]
t,e = dfe.iloc[:,0], dfe.iloc[:,1]

Inp = input("Enter s for solution plot \n or \n Enter e for error plot\n")


if(Inp=="s"):
    plt.plot(t,y,label="Solution")
else:
    plt.plot(t,e,label="error")

plt.legend()
plt.show()
