"""L2-genome gated synteny repair (GSR): unit tests of the rule and of its default-off key."""
import numpy as np

from ancestra.merge import MergeConfig, _gsr_repair


def _run(cfg, sup, cs=None):
    # units: 0,1 (child 0) grouped with 2 (child 1); 3 (child 1) grouped with 4 (child 0);
    # unit 5 (child 1) is alone. Candidate edges from 5: to 0 (group A) and to 4 (group B).
    child = np.array([0, 0, 1, 1, 0, 1])
    labels = np.array([0, 0, 0, 3, 3, 5])
    ca = np.array([0, 4, 0])
    cb = np.array([5, 5, 2])
    cs = np.array([0.8, 0.9, 0.95]) if cs is None else cs
    st = {}
    out = _gsr_repair(6, child, ca, cb, cs, np.asarray(sup, float), labels, cfg, st)
    return out, st


def test_default_off_not_in_key():
    k = MergeConfig().key()
    assert not any(x.startswith("gsr") for x in k)


def test_joins_unique_best_supported_group_lacking_child():
    cfg = MergeConfig(gsr="any")
    # group A (label 0) already holds child 1 (unit 2) -> not allowed; group B holds child 1 (unit 3)
    out, st = _run(cfg, [3, 0, 0])
    assert st["gsr_repaired"] == 0 and out[5] == 5


def test_join_when_group_lacks_child():
    cfg = MergeConfig(gsr="any")
    child = np.array([0, 0, 1])
    labels = np.array([0, 0, 2])
    ca, cb = np.array([0, 1]), np.array([2, 2])
    st = {}
    out = _gsr_repair(3, child, ca, cb, np.array([0.9, 0.8]), np.array([3.0, 0.0]), labels, cfg, st)
    assert out[2] == 0 and st["gsr_repaired"] == 1


def test_margin_and_sequence_gate():
    cfg = MergeConfig(gsr="any")
    child = np.array([0, 1, 2])
    labels = np.array([0, 1, 2])      # all alone: unit 2 may join 0 or 1
    ca, cb = np.array([0, 1]), np.array([2, 2])
    st = {}
    # tie in support -> margin fails
    out = _gsr_repair(3, child, ca, cb, np.array([0.9, 0.9]), np.array([2.0, 2.0]), labels, cfg, st, np.array([0, 0, 1], bool))
    assert st["gsr_repaired"] == 0
    # support points to a partner below rho x best -> sequence gate fails
    st = {}
    out = _gsr_repair(3, child, ca, cb, np.array([0.5, 0.9]), np.array([3.0, 0.0]), labels, cfg, st, np.array([0, 0, 1], bool))
    assert st["gsr_repaired"] == 0
    st = {}
    out = _gsr_repair(3, child, ca, cb, np.array([0.8, 0.9]), np.array([3.0, 0.0]), labels, cfg, st,
                      np.array([0, 0, 1], bool))
    assert st["gsr_repaired"] == 1 and out[2] == out[0]
