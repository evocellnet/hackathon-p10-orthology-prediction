"""Leakage control (DESIGN.md §4): method code may read only declared inputs.

The method package ``ancestra/`` must not reference evaluation data: QfO reference
data / gold standards, OrthoBench RefOGs, SwissTree, TreeFam, the QfO lineage tree and
species-tree samples, published participant predictions, VGNC, or the benchmark
harness (``bench_repo``). Such references are allowed only in ``docs/``, ``tests/`` and
``scripts/`` (evaluation code).

1. ``test_static_scan``: every Python file under ``ancestra/`` is tokenised; identifiers,
   attribute names, import names and string literals (docstrings excluded) are matched
   case-insensitively against FORBIDDEN. Comments and docstrings are provenance notes
   ("selected on OrthoBench dev RefOGs") and not data access; they are listed by
   ``test_static_scan_report`` but do not fail. Non-Python files in ``ancestra/`` (schemas)
   are scanned in full.
2. ``test_runtime_audit``: the pipeline (default S0 families -> scores -> merge -> output)
   runs on a toy simulation under a ``sys.addaudithook`` hook that records every file
   opened by Python and every subprocess command line; none may contain FORBIDDEN.
"""
import io
import os
import re
import sys
import tokenize

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKG = os.path.join(ROOT, "ancestra")
FORBIDDEN = ("reference_data", "refog", "swisstree", "treefam", "lineage_tree",
             "species_tree_samples", "participants", "vgnc", "bench_repo")
PAT = re.compile("|".join(re.escape(f) for f in FORBIDDEN), re.IGNORECASE)


def _py_files():
    for dp, dn, fn in os.walk(PKG):
        dn[:] = [d for d in dn if d != "__pycache__"]
        for f in sorted(fn):
            yield os.path.join(dp, f)


def _docstring_lines(src: str) -> set:
    """Line numbers of module/class/function docstrings."""
    import ast
    out = set()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                out.update(range(body[0].lineno, body[0].end_lineno + 1))
    return out


def scan_file(path: str):
    """Return (code_hits, note_hits): lists of (path:line, token text)."""
    rel = os.path.relpath(path, ROOT)
    with open(path, "rb") as fh:
        raw = fh.read()
    if not path.endswith(".py"):
        try:
            text = raw.decode()
        except UnicodeDecodeError:
            return [], []
        hits = [(f"{rel}:{i}", m.group(0)) for i, line in enumerate(text.splitlines(), 1)
                for m in PAT.finditer(line)]
        return hits, []
    src = raw.decode()
    doc = _docstring_lines(src)
    code, notes = [], []
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type not in (tokenize.NAME, tokenize.STRING, tokenize.COMMENT, tokenize.OP) \
                and not (sys.version_info >= (3, 12) and tok.type == getattr(tokenize, "FSTRING_MIDDLE", -1)):
            continue
        if not PAT.search(tok.string):
            continue
        where = f"{rel}:{tok.start[0]}"
        if tok.type == tokenize.COMMENT or (tok.type == tokenize.STRING and tok.start[0] in doc):
            notes.append((where, tok.string.strip()[:80]))
        else:
            code.append((where, tok.string[:80]))
    return code, notes


def test_forbidden_pattern_self_check():
    """The scanner must catch the forbidden tokens (positive control) and ignore clean code."""
    for bad in ("x = 'data/qfo/reference_data/2022.1/mapping.json.gz'", "from foo import read_refogs",
                "p = os.path.join(d, 'RefOG001.txt')", "t = load('SwissTree')", "v = VGNC_FILE",
                "f = 'participants_2022'", "s = 'species_tree_samples_Euk.nwk'"):
        assert PAT.search(bad), bad
    for ok in ("x = load_proteomes(d)", "tree = load_guide_tree(p, sp)", "fam = families.npy"):
        assert not PAT.search(ok), ok


def test_static_scan():
    hits = []
    for f in _py_files():
        code, _ = scan_file(f)
        hits += code
    msg = "\n".join(f"  {w}: {t}" for w, t in hits)
    assert not hits, ("forbidden evaluation references in method code (ancestra/); move evaluation "
                      "helpers to scripts/ or tests/:\n" + msg)


def test_static_scan_report(capsys):
    """Comment/docstring mentions: informational only (provenance), printed with -s."""
    notes = []
    for f in _py_files():
        notes += scan_file(f)[1]
    with capsys.disabled():
        if notes:
            print(f"\n[test_no_leak] {len(notes)} comment/docstring mentions (allowed):")
            for w, t in notes[:20]:
                print(f"  {w}: {t}")


def test_runtime_audit(sim_basic, tmp_path):
    from ancestra.pipeline import run_pipeline
    _, simdir = sim_basic
    seen = []
    state = {"on": True}

    def hook(event, args):
        if not state["on"]:
            return
        pathlike = bool(args) and isinstance(args[0], (str, bytes, os.PathLike))
        if event == "open" and pathlike:
            seen.append(("open", os.fsdecode(args[0])))
        elif event == "subprocess.Popen" and len(args) > 1:
            seen.append(("exec", " ".join(map(str, args[1] or []))))
        elif event in ("os.listdir", "os.scandir") and pathlike:
            seen.append(("list", os.fsdecode(args[0])))

    sys.addaudithook(hook)
    try:
        run_pipeline(os.path.join(simdir, "proteomes"), os.path.join(simdir, "species_tree.nwk"),
                     str(tmp_path / "run"), threads=4, log=lambda *_: None)
    finally:
        state["on"] = False          # audit hooks cannot be removed; disable it
    assert any(k == "exec" and "diamond" in c for k, c in seen), "audit hook saw no tool calls"
    bad = [(k, p) for k, p in seen if PAT.search(p)]
    assert not bad, f"pipeline touched evaluation data: {bad[:10]}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-s"]))
