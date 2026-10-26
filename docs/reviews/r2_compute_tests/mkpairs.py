import numpy as np, random
from lxml import etree
O='{http://orthoXML.org/2011/}'
f='/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/data/qfo/participants_2022/FastOMA_v0.1.4_qfo2022_17Nov3.orthoxml'
id2acc={}
fams=[]
for ev,el in etree.iterparse(f,events=('end',),tag=(O+'gene',O+'groups')):
    if el.tag==O+'gene':
        id2acc[el.get('id')]=el.get('protId'); 
    else:
        for top in el:
            if top.tag in (O+'orthologGroup',O+'paralogGroup'):
                fams.append([g.get('id') for g in top.iter(O+'geneRef')])
acc2key={}
for l in open('qdb.lookup'):
    k,a,_=l.split('\t'); acc2key[a]=int(k)
random.seed(11); random.shuffle(fams)
tot=0; out=[]; famfa=[]
for fam in fams:
    keys=[acc2key.get(id2acc[g].split('|')[-1] if '|' in id2acc[g] else id2acc[g]) for g in fam]
    keys=[k for k in keys if k is not None]
    if len(keys)<2 or len(keys)>1500: continue
    out.append(keys); tot+=len(keys)**2
    if tot>5e6: break
print('fams',len(out),'genes',sum(map(len,out)),'pairs',tot, 'miss check sample', list(id2acc.values())[:2])
with open('pairs.tsv','w') as fh:
    for keys in out:
        for q in keys:
            for t in keys:
                fh.write(f'{q}\t{t}\t0\t0\n')
with open('subset.keys','w') as fh:
    for keys in out:
        for q in keys: fh.write(f'{q}\n')
