import numpy as np, pandas as pd
from lxml import etree
from collections import Counter, defaultdict
dc=pd.read_csv('dc.tsv',sep='\t',header=None,names=['c','m'])
acc=lambda s: s.split('|')[1] if '|' in s else s
dc['c']=dc.c.map(acc); dc['m']=dc.m.map(acc)
s=dc.groupby('c').size().sort_values(ascending=False).values.astype(float)
N=len(dc)
print('clusters',len(s),'singletons',(s==1).sum(),'genes in size>=2',int(s[s>=2].sum()),'top10',s[:10].astype(int).tolist())
print('sum g^2',f'{(s**2).sum():.3e}','N^2',f'{N*N:.3e}','ratio',f'{N*N/(s**2).sum():.0f}x')
for cap in (1000,2000,5000): print('cap',cap,f'{np.minimum(s**2,cap*s).sum():.3e}','n>cap',(s>cap).sum())
# fission vs FastOMA root HOGs
m2c=dict(zip(dc.m,dc.c))
O='{http://orthoXML.org/2011/}'
f='/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/data/qfo/participants_2022/FastOMA_v0.1.4_qfo2022_17Nov3.orthoxml'
id2acc={}; fams=[]
for ev,el in etree.iterparse(f,events=('end',),tag=(O+'gene',O+'groups')):
    if el.tag==O+'gene': id2acc[el.get('id')]=acc(el.get('protId'))
    else:
        for top in el:
            if top.tag in (O+'orthologGroup',O+'paralogGroup'):
                fams.append([id2acc[g.get('id')] for g in top.iter(O+'geneRef')])
whole=0; frac=[]; nclus=[]
for fam in fams:
    cs=Counter(m2c.get(a) for a in fam)
    top=cs.most_common(1)[0][1]
    frac.append(top/len(fam)); nclus.append(len(cs)); whole+= (len(cs)==1)
frac=np.array(frac); nclus=np.array(nclus); sz=np.array([len(f) for f in fams])
print('FastOMA rootHOGs',len(fams),'fully in one DC cluster',f'{whole/len(fams):.3f}','gene-weighted mean largest-share',f'{(frac*sz).sum()/sz.sum():.3f}','median #clusters per HOG',np.median(nclus),'HOGs>=10 genes: frac whole',f'{((nclus==1)&(sz>=10)).sum()/(sz>=10).sum():.3f}')
# fusion: how many FastOMA HOGs per DC cluster (top clusters)
h2=defaultdict(set)
for i,fam in enumerate(fams):
    for a in fam: h2[m2c.get(a)].add(i)
nh=np.array([len(v) for v in h2.values()])
print('DC clusters spanning >1 FastOMA HOG',(nh>1).sum(),'max HOGs in one cluster',nh.max())
