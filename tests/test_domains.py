"""L2-domains: Pfam architecture parsing, shuffled control, unit architectures, merge hooks (default off)."""
import numpy as np

from ancestra.domains import DomainArch, parse_domtbl
from ancestra.merge import MergeConfig, NodeContext, decide_node, run_merge, to_orthoxml
from ancestra.orthoxml import hogs_to_pairs


def _domline(prot, acc, ie, s, e):
    f = [prot, "-", "500", "dom", acc + ".1", "100", "1e-10", "50", "0", "1", "1", "1e-10", str(ie), "50", "0",
         "1", "100", "1", "1", str(s), str(e), "0.9", "desc"]
    return " ".join(f)


def test_parse_overlap_clan_collapse(tmp_path):
    p = tmp_path / "x.domtbl"
    p.write_text("\n".join([
        "# comment",
        _domline("P1", "PF00001", 1e-20, 10, 100),
        _domline("P1", "PF00002", 1e-5, 20, 110),    # overlaps the better hit > 50% -> dropped
        _domline("P1", "PF00400", 1e-8, 150, 190),   # WD40 repeat x3, collapsed
        _domline("P1", "PF00400", 1e-8, 200, 240),
        _domline("P1", "PF00400", 1e-8, 250, 290),
        _domline("P2", "PF00003", 1e-8, 5, 50),
    ]) + "\n")
    a = parse_domtbl(str(p), {"PF00400": "CL0186"})
    assert a["P1"] == ("PF00001", "CL0186")
    assert a["P2"] == ("PF00003",)


def test_shuffle_keeps_species_and_annotation():
    names = [f"g{i}" for i in range(40)]
    sp = np.array([i % 2 for i in range(40)])
    arch = {f"g{i}": (("A",) if i % 4 == 0 else ("B", "C")) for i in range(40) if i % 5}
    real = DomainArch(names, sp, arch)
    sh = DomainArch(names, sp, arch, shuffle_seed=11)
    assert ((real.code >= 0) == (sh.code >= 0)).all()
    for s in (0, 1):
        m = (sp == s) & (real.code >= 0)
        assert sorted(real.archs[c] for c in real.code[m]) == sorted(sh.archs[c] for c in sh.code[m])
    assert (real.code != sh.code).any()


def test_unit_codes_mode_and_jaccard():
    names = ["a", "b", "c", "d", "e"]
    arch = {"a": ("X", "Y"), "b": ("X", "Y"), "c": ("X",), "d": ("Z",)}
    da = DomainArch(names, np.zeros(5, int), arch)
    g_sorted = np.array([0, 1, 2, 3, 4])
    b = np.array([0, 3, 4, 5])              # unit0 = a,b,c ; unit1 = d ; unit2 = e (unannotated)
    ctx = NodeContext(node="x", children=None, units=None, child=np.array([0, 1, 1]),
                      cand_a=np.array([0, 0]), cand_b=np.array([1, 2]), cand_score=np.ones(2), anchors=None,
                      adjacency=None, members=lambda i, _g=g_sorted, _b=b: _g[_b[i]:_b[i + 1]])
    uc = da.unit_codes(ctx)
    assert da.archs[uc[0]] == ("X", "Y") and da.archs[uc[1]] == ("Z",) and uc[2] == -1
    j, same = da.edge_sim(ctx)
    assert j[0] == 0.0 and np.isnan(j[1]) and not same.any()


class _Arch:
    """Fake arch object: fixed per-edge answers."""
    def __init__(self, j, same):
        self.j, self.same = np.asarray(j, float), np.asarray(same, bool)

    def edge_sim(self, ctx, sel=None):
        sel = np.arange(len(self.j)) if sel is None else sel
        return self.j[sel], self.same[sel]


def test_gate_and_veto_in_decide_node():
    # child0: a ; child1: b1, b2. a-b1 0.50, a-b2 0.42 (0.84 of best: mutual best only at tol 0.8)
    child = np.array([0, 1, 1])
    lu, lv, sc = np.array([0, 0]), np.array([1, 2]), np.array([0.50, 0.42])
    base = MergeConfig(tree_mode=False)
    info = {"members": lambda i: np.array([i])}
    _, s0 = decide_node(3, child, lu, lv, sc, base, node_info=info, arch=_Arch([1, 1], [True, True]))
    assert s0["n_mb"] == 1 and "arch_gate_added" not in s0          # mode off: arch ignored
    _, s1 = decide_node(3, child, lu, lv, sc, MergeConfig(tree_mode=False, arch_mode="gate"), node_info=info,
                        arch=_Arch([1, 1], [True, False]))
    assert s1.get("arch_gate_added", 0) == 0                         # different architecture: no rescue
    _, s2 = decide_node(3, child, lu, lv, sc, MergeConfig(tree_mode=False, arch_mode="gate"), node_info=info,
                        arch=_Arch([1, 1], [True, True]))
    assert s2["arch_gate_added"] == 1
    _, s3 = decide_node(3, child, lu, lv, sc, MergeConfig(tree_mode=False, arch_mode="veto"), node_info=info,
                        arch=_Arch([0.0, 1], [False, True]))
    assert s3["arch_veto"] == 1
    assert "arch_mode" not in MergeConfig().key()                    # default stage hash unchanged


def test_run_merge_arch_off_identical():
    from test_merge import GENES, SCORES, TREE, make_case
    gt, st, tree = make_case(GENES, SCORES, TREE)
    arch = DomainArch(gt.names, gt.species_of, {g: ("X",) for g in gt.names}, shuffle_seed=11)
    r0 = run_merge(gt, tree, st, MergeConfig())
    r1 = run_merge(gt, tree, st, MergeConfig(), arch=arch)
    assert hogs_to_pairs(to_orthoxml(r0, gt, tree)) == hogs_to_pairs(to_orthoxml(r1, gt, tree))
    r2 = run_merge(gt, tree, st, MergeConfig(arch_mode="gate"), arch=arch)   # runs end to end
    assert len(hogs_to_pairs(to_orthoxml(r2, gt, tree))) >= len(hogs_to_pairs(to_orthoxml(r0, gt, tree)))
