import gzip,sys,pandas as pd,re
sp={'Homo_sapiens':'Homo_sapiens.GRCh38.99.gff3.gz','Pan_troglodytes':'Pan_troglodytes.Pan_tro_3.0.99.gff3.gz','Mus_musculus':'Mus_musculus.GRCm38.99.gff3.gz','Rattus_norvegicus':'Rattus_norvegicus.Rnor_6.0.99.gff3.gz','Canis_familiaris':'Canis_familiaris.CanFam3.1.99.gff3.gz','Monodelphis_domestica':'Monodelphis_domestica.ASM229v1.99.gff3.gz'}
for s,f in sp.items():
    rows=[]
    for line in gzip.open('data/orthobench/ensembl_99/gff3/'+f,'rt'):
        if line[0]=='#': continue
        x=line.split('\t')
        if len(x)<9 or x[2] not in ('gene','ncRNA_gene','pseudogene'): continue
        a=dict(kv.split('=',1) for kv in x[8].strip().split(';') if '=' in kv)
        rows.append((x[0],int(x[3]),int(x[4]),x[2],a.get('biotype',''),a.get('Name',''),a.get('gene_id',a.get('ID','').replace('gene:',''))))
    d=pd.DataFrame(rows,columns=['contig','start','end','type','biotype','name','gene_id'])
    d.to_parquet(f'work/S8/gff_genes_{s}.parquet'); print(s,len(d),d.biotype.value_counts().head(5).to_dict())
