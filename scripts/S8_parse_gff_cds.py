# CDS segments of the input-proteome proteins (target-side overlap test for S8)
import gzip,pandas as pd
sp={'Homo_sapiens':'Homo_sapiens.GRCh38.99.gff3.gz','Pan_troglodytes':'Pan_troglodytes.Pan_tro_3.0.99.gff3.gz','Mus_musculus':'Mus_musculus.GRCm38.99.gff3.gz','Rattus_norvegicus':'Rattus_norvegicus.Rnor_6.0.99.gff3.gz','Canis_familiaris':'Canis_familiaris.CanFam3.1.99.gff3.gz','Monodelphis_domestica':'Monodelphis_domestica.ASM229v1.99.gff3.gz'}
for s,f in sp.items():
    want=set(pd.read_parquet(f'data/genomes/orthobench/{s}.parquet',columns=['protein_id']).protein_id)
    rows=[]
    for line in gzip.open('data/orthobench/ensembl_99/gff3/'+f,'rt'):
        if line[0]=='#': continue
        x=line.split('\t',8)
        if len(x)<9 or x[2]!='CDS': continue
        i=x[8].find('protein_id=')
        if i<0: continue
        p=x[8][i+11:].split(';')[0].strip()
        if p in want: rows.append((p,x[0],int(x[3]),int(x[4]),x[6]))
    d=pd.DataFrame(rows,columns=['protein_id','contig','start','end','strand'])
    d.to_parquet(f'work/S8/cds_{s}.parquet'); print(s,len(d),d.protein_id.nunique(),len(want))
