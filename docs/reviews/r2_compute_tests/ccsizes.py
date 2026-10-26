import sys, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
df=pd.read_csv(sys.argv[1],sep='\t',header=None,usecols=[0,1],names=['a','b'],dtype=str)
ids=pd.Index(pd.unique(np.concatenate([df.a.values,df.b.values])))
a=ids.get_indexer(df.a); b=ids.get_indexer(df.b); n=len(ids)
M=sp.coo_matrix((np.ones(len(a),dtype=np.int8),(a,b)),shape=(n,n))
k,lab=connected_components(M,directed=False)
s=np.sort(np.bincount(lab))[::-1].astype(float)
print('genes',n,'comps',k,'top10',s[:10].astype(int).tolist())
print('sum g^2',f'{(s**2).sum():.3e}','top1 share',f'{s[0]**2/(s**2).sum():.2f}')
for cap in (1000,2000,5000):
    print('cap',cap,f'{np.minimum(s**2,cap*s).sum():.3e}','n>cap',(s>cap).sum(),'genes>cap',int(s[s>cap].sum()))
