#!/usr/bin/env python
"""Merge-only variant of a finished `python -m ancestra run`, through the unmodified frozen CLI.

    python run_cached.py BASE_RUN -- <`ancestra run` args (without 'run') incl. --out NEW_RUN>

Before the CLI runs, NEW_RUN/stages/{input,families,scores} are hard-linked from BASE_RUN (as in
scripts/rescue_ob_arms.sh). For a *different guide tree* (`--tree`) the input-stage key changes (it
contains the tree md5) and therefore the keys of every downstream stage. The families (S0) and
scores (S2) stages do not read the guide tree (pipeline.run_pipeline: s_fam uses only the gene table
and FamilyConfig; s_scores = score_with_rescue(gt, fam, scoring)), so for those two stages this
wrapper accepts the base run's checkpoint under the new key: when the CLI asks for stage `families`
or `scores` and the linked DONE holds a different key, the DONE file is replaced (unlinked first,
so the base run's file is untouched) by the new key, and the stage is then reported as `skipped`.
`input`, `rescue`, `merge` and `output` always run with their own keys (input is rerun -> writes the
new guide tree). Nothing in ancestra/ is modified; the wrapper only patches `Runner.stage` in-process.
Every adoption is logged in NEW_RUN/reused_stages.json.
"""
import json
import os
import shutil
import subprocess
import sys

ROOT = "/scratch/spascare/agentic-ai/workspaces/p10-orthology-prediction"
sys.path.insert(0, ROOT)

import ancestra.pipeline as P  # noqa: E402

TREE_FREE = ("families", "scores")


def main():
    i = sys.argv.index("--")
    base = os.path.abspath(sys.argv[1])
    args = sys.argv[i + 1:]
    out = args[args.index("--out") + 1]
    os.makedirs(os.path.join(out, "stages"), exist_ok=True)
    for s in ("input",) + TREE_FREE:
        dst = os.path.join(out, "stages", s)
        if not os.path.exists(dst):
            subprocess.run(["cp", "-al", os.path.join(base, "stages", s), dst], check=True)
    log = []
    orig = P.Runner.stage

    def stage(self, name, key, fn):
        if name in TREE_FREE:
            d = os.path.join(self.sdir, name)
            done = os.path.join(d, "DONE")
            h = P._hash(key)
            bdone = os.path.join(base, "stages", name, "DONE")
            if os.path.exists(done) and open(done).read().strip() != h:
                old = open(done).read().strip()
                assert old == open(bdone).read().strip(), "stage dir is not the base run's checkpoint"
                os.unlink(done)                      # break the hard link; base DONE unchanged
                with open(done, "w") as fh:
                    fh.write(h + "\n")
                log.append(dict(stage=name, base_key=old, new_key=h, base=base))
        return orig(self, name, key, fn)

    P.Runner.stage = stage
    sys.argv = ["run_cached.py(base=%s)" % base] + args
    from ancestra.__main__ import main as amain
    try:
        rc = amain(["run"] + args)
    finally:
        json.dump(dict(base=base, argv=args, adopted=log), open(os.path.join(out, "reused_stages.json"), "w"),
                  indent=1)
    sys.exit(rc)


if __name__ == "__main__":
    main()
