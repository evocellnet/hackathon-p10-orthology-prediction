#!/usr/bin/env bash
# Smoke test of all tools in conda_envs/orth. Usage: bash scripts/smoke_test_tools.sh
set -u
ROOT=/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction
source /opt/conda/etc/profile.d/conda.sh && conda activate $ROOT/conda_envs/orth
export HF_HOME=$ROOT/data/models HF_HUB_OFFLINE=1
T=$(mktemp -d); cd $T; fail=0
cat > a.fa <<F
>p1
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV
>p2
MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVA
>p3
MKSAYIAKQRQISFVKSHFSRQLEERLGLIEVQPPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV
>p4
MNTAYLAKQRQISFVKSHFSRQLEQRLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDGQFEVV
F
run(){ n=$1; shift; if "$@" >log.$n 2>&1; then echo "PASS $n"; else echo "FAIL $n"; tail -3 log.$n; fail=1; fi; }
run diamond bash -c 'diamond makedb --in a.fa -d a -p 2 && diamond blastp -q a.fa -d a -o d.tsv -p 2 && test -s d.tsv'
run mmseqs bash -c 'mmseqs easy-search a.fa a.fa m.tsv tmp --threads 2 && test -s m.tsv'
run foldseek bash -c 'foldseek version'
run mafft bash -c 'mafft --quiet --auto a.fa > a.aln && test -s a.aln'
run famsa bash -c 'famsa -t 2 a.fa f.aln && test -s f.aln'
run muscle bash -c 'muscle -align a.fa -output mu.aln && test -s mu.aln'
run fasttree bash -c 'FastTree -quiet a.aln > ft.nwk && test -s ft.nwk'
run iqtree bash -c 'iqtree -s a.aln -m LG -T 1 -pre iq -quiet && test -s iq.treefile'
run trimal bash -c 'trimal -in a.aln -out t.aln -automated1 && test -s t.aln'
run hmmer bash -c 'hmmbuild h.hmm a.aln && hmmsearch --tblout h.tbl h.hmm a.fa && test -s h.tbl'
run blast bash -c 'makeblastdb -in a.fa -dbtype prot -out bdb && blastp -query a.fa -db bdb -out b.tsv -outfmt 6 && test -s b.tsv'
run cd-hit bash -c 'cd-hit -i a.fa -o c.fa -c 0.9 && test -s c.fa'
run mcl bash -c 'printf "a b 1\nb c 1\nc a 1\nd e 1\n" > g.abc && mcl g.abc --abc -o g.mcl && test -s g.mcl'
run orthofinder bash -c 'orthofinder -h'
run omamer bash -c 'omamer -h'
run fastoma bash -c 'fastoma-core -h || fastoma-helper -h || fastoma-batch-fasta -h; ls $CONDA_PREFIX/bin | grep -qi fastoma'
run sonicparanoid bash -c 'sonicparanoid -h'
run proteinortho bash -c 'proteinortho -v'
run generax bash -c 'generax --help'
run aster bash -c 'astral-pro3 -h || astral4 -h || ls $CONDA_PREFIX/bin | grep -iE "astral|wastral"'
run newick_utils bash -c 'echo "((a:1,b:1):1,c:2);" | nw_topology - | grep -q a'
run snakemake bash -c 'snakemake --version'
run pyimports python -c "
import numpy,scipy,pandas,numba,sklearn,networkx,ete3,Bio,dendropy,pyham,lxml,faiss,pyarrow,matplotlib,seaborn,tqdm,pytest,igraph,leidenalg,torch,transformers,sentencepiece,esm
print('ok')"
run esm2 python -c "
from transformers import AutoTokenizer, AutoModel
for r in ['facebook/esm2_t6_8M_UR50D','facebook/esm2_t12_35M_UR50D']:
    AutoModel.from_pretrained(r)(**AutoTokenizer.from_pretrained(r)('MKTAYIAKQR',return_tensors='pt'))"
run prostt5tok python -c "
from transformers import AutoTokenizer
AutoTokenizer.from_pretrained('Rostlab/ProstT5',do_lower_case=False)('<AA2fold> A E T')"
echo; [ $fail = 0 ] && echo "ALL PASS" || echo "SOME FAILED"; rm -rf $T; exit $fail
