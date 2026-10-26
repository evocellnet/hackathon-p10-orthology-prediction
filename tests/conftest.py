"""Shared fixtures: small simulations written once per test session."""
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ancestra.simulate import SimParams, simulate, write_simulation  # noqa: E402

QFO_DIR = os.path.join(ROOT, "data", "qfo", "benchmark-webservice")
ORTHOBENCH_DIR = os.path.join(ROOT, "data", "orthobench", "Open_Orthobench", "BENCHMARKS")
ORTHOBENCH_POSITIVE = os.path.join(ROOT, "data", "orthobench", "controls", "positive_refogs.txt")


@pytest.fixture(scope="session")
def sim_basic(tmp_path_factory):
    """8 species x 80 genes, default model (no HGT, no WGD)."""
    res = simulate(SimParams(n_species=8, n_genes=80), seed=11)
    out = str(tmp_path_factory.mktemp("sim_basic"))
    write_simulation(res, out)
    return res, out


@pytest.fixture(scope="session")
def sim_complex(tmp_path_factory):
    """6 species x 60 genes with HGT, a WGD and neofunctionalisation."""
    p = SimParams(n_species=6, n_genes=60, hgt_rate=0.08, wgd_branch="random", neo_factor=3.0,
                  dup_rate=0.3, loss_rate=0.3)
    res = simulate(p, seed=23)
    out = str(tmp_path_factory.mktemp("sim_complex"))
    write_simulation(res, out)
    return res, out
