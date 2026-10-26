# Estimate ancestral proteome sizes per node of QfO lineage tree using FastOMA HOGs as a proxy
import sys
from lxml import etree
from collections import defaultdict
Q='/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction/data/qfo/'
ns={'p':'http://www.phyloxml.org'}
t=etree.parse(Q+'reference_data/2022.1/lineage_tree.phyloxml')
root=t.find('.//p:phylogeny/p:clade',ns)
nodes=[]  # (name, leafset, child leafsets)
def walk(c):
    kids=c.findall('p:clade',ns)
    if not kids:
        return frozenset([c.find('p:taxonomy/p:code',ns).text])
    ks=[walk(k) for k in kids]
    L=frozenset().union(*ks)
    nodes.append((c.find('p:taxonomy/p:scientific_name',ns).text,L,ks))
    return L
walk(root)
O='{http://orthoXML.org/2011/}'
gene2sp={}; spsize=defaultdict(int)
# pass 1: genes
for ev,el in etree.iterparse(Q+'participants_2022/FastOMA_v0.1.4_qfo2022_17Nov3.orthoxml',events=('end',),tag=O+'species'):
    sp=el.get('name')
    for g in el.iter(O+'gene'):
        gene2sp[g.get('id')]=sp; spsize[sp]+=1
    el.clear()
print('genes',len(gene2sp),'species',len(spsize),file=sys.stderr)
# pass 2: for each top-level group, compute nested structure: list of (species set, member genes) for orthologGroups
# For node v, units within a top group = number of maximal orthologGroups (or genes) whose species set subset of L(v), restricted to genes of L(v)
anc=defaultdict(int)
ingroup=set()
def proc(el):
    # returns (speciesset, count dict per node) - recursion: units(v) for subtree element
    tag=el.tag.replace(O,'')
    if tag=='geneRef':
        gid=el.get('id'); ingroup.add(gid); sp=gene2sp[gid]
        return frozenset([sp]), {i:(1 if sp in L else 0) for i,(n,L,_) in enumerate(nodes)}
    ch=[proc(c) for c in el if c.tag.replace(O,'') in ('geneRef','orthologGroup','paralogGroup')]
    S=frozenset().union(*[c[0] for c in ch]) if ch else frozenset()
    cnt={}
    for i,(n,L,_) in enumerate(nodes):
        if tag=='orthologGroup' and S and S<=L: cnt[i]=1
        else: cnt[i]=sum(c[1][i] for c in ch)
    return S,cnt
ctx=etree.iterparse(Q+'participants_2022/FastOMA_v0.1.4_qfo2022_17Nov3.orthoxml',events=('end',),tag=O+'groups')
for ev,gr in ctx:
    for top in gr:
        if top.tag.replace(O,'') not in ('orthologGroup','paralogGroup'): continue
        S,cnt=proc(top)
        for i,v in cnt.items(): anc[i]+=v
# singletons
for gid,sp in gene2sp.items():
    if gid not in ingroup:
        for i,(n,L,_) in enumerate(nodes):
            if sp in L: anc[i]+=1
def size(L):
    if len(L)==1: return spsize[next(iter(L))]
    for i,(n,LL,_) in enumerate(nodes):
        if LL==L: return anc[i]
tot=0;tot_pairs=0
print('node\tnleaves\tnchild\tanc_size\tsum_leafgenes\tchild_sizes\tcross_pairs')
for i,(n,L,ks) in enumerate(nodes):
    cs=[size(k) for k in ks]
    cross=sum(cs[a]*cs[b] for a in range(len(cs)) for b in range(a+1,len(cs)))
    tot_pairs+=cross
    print(f"{n}\t{len(L)}\t{len(ks)}\t{anc[i]}\t{sum(spsize[s] for s in L)}\t{cs}\t{cross:.3e}")
N=sum(spsize.values())
print('TOTAL cross pairs (k=1)',f'{tot_pairs:.3e}','all-vs-all N^2/2',f'{N*N/2:.3e}', 'N',N)
