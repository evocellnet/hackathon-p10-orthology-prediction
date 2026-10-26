"""S0 determinism: identical families regardless of --threads, input order and ID names."""
import os
import random

import pytest

from ancestra.families import FamilyConfig, read_families_tsv, run_families
from ancestra.pipeline import DEFAULT_FAMILY_CONFIG


def _read(fa):
    recs, name = [], None
    for line in open(fa):
        if line.startswith(">"):
            name = line[1:].split()[0]
            recs.append([name, ""])
        else:
            recs[-1][1] += line.strip()
    return recs


def _partition(m, rename=None):
    groups = {}
    for g, f in m.items():
        groups.setdefault(f, set()).add(rename(g) if rename else g)
    return {frozenset(x) for x in groups.values()}


@pytest.mark.parametrize("gmax", [None, 8])           # 8: forces Leiden splits + reassignment
def test_s0_identical_across_threads_order_and_names(sim_basic, tmp_path, gmax):
    _, simdir = sim_basic
    recs = []
    pdir = os.path.join(simdir, "proteomes")
    for f in sorted(os.listdir(pdir)):
        recs += _read(os.path.join(pdir, f))
    fa1 = tmp_path / "in1.fa"
    fa1.write_text("".join(f">{n}\n{s}\n" for n, s in recs))
    shuf = list(recs)
    random.Random(5).shuffle(shuf)
    fa2 = tmp_path / "in2.fa"                          # shuffled order, renamed IDs, "*"
    fa2.write_text("".join(f">X_{n}\n{s}*\n" for n, s in shuf))   # + terminal stop codons
    parts = []
    for fa, t, tag in ((fa1, 2, "a"), (fa1, 8, "b"), (fa2, 8, "c")):
        cfg = FamilyConfig(**DEFAULT_FAMILY_CONFIG)
        cfg.threads = t
        cfg.s2_gmax = gmax
        out = tmp_path / f"run_{tag}"
        st = run_families(str(fa), str(out), cfg)
        m = read_families_tsv(str(out / "families.tsv"))
        parts.append(_partition(m, (lambda g: g[2:]) if tag == "c" else None))
        if gmax:
            assert st["n_components_split"] >= 1
    assert parts[0] == parts[1] == parts[2]
    # family IDs (labels) are identical too for the same input
    a = read_families_tsv(str(tmp_path / "run_a" / "families.tsv"))
    b = read_families_tsv(str(tmp_path / "run_b" / "families.tsv"))
    assert a == b
