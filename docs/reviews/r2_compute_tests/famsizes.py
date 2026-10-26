import sys, numpy as np
from lxml import etree
f=sys.argv[1]; O='{http://orthoXML.org/2011/}'
sizes=[]
for ev,gr in etree.iterparse(f,events=('end',),tag=O+'groups'):
    for top in gr:
        if top.tag in (O+'orthologGroup',O+'paralogGroup'):
            sizes.append(sum(1 for _ in top.iter(O+'geneRef')))
s=np.array(sorted(sizes))[::-1]
N=988778
print('groups',len(s),'genes_in',s.sum(),'singletons',N-s.sum())
print('top10',s[:10].tolist())
print('sum g^2 (ordered pairs)',f'{(s.astype(float)**2).sum():.3e}','share top1/top10/top100', *(f'{(s[:k].astype(float)**2).sum()/(s.astype(float)**2).sum():.2f}' for k in (1,10,100)))
for cap in (500,1000,2000,5000):
    print('cap',cap,'sum min(g^2,cap*g)',f'{np.minimum(s.astype(float)**2,cap*s.astype(float)).sum():.3e}','n>cap',(s>cap).sum(),'genes>cap',s[s>cap].sum())
