#!/bin/bash
# Key S3/S5 arms on one OrthoBench run directory (dev only): off / constrained real / constrained
# shuffled (3 seeds) / locus repair (real + shuffled). Usage: ob_key_arms.sh RUN_DIR TAG
R=$1; T=$2
C="--syn-constrain --set ctx_veto=true ctx_r3=true"
python scripts/ob_merge_iter.py $R --arm default --label "$T off" &
python scripts/ob_merge_iter.py $R --arm default --synteny real $C --label "$T syn real constrained+veto+r3" &
python scripts/ob_merge_iter.py $R --arm default --synteny shuffled --syn shuffle_seed=11 $C --label "$T syn shuffled(s11) constrained+veto+r3" &
wait
python scripts/ob_merge_iter.py $R --arm default --synteny shuffled --syn shuffle_seed=23 $C --label "$T syn shuffled(s23) constrained+veto+r3" &
python scripts/ob_merge_iter.py $R --arm default --synteny shuffled --syn shuffle_seed=47 $C --label "$T syn shuffled(s47) constrained+veto+r3" &
python scripts/ob_merge_iter.py $R --arm default --synteny real $C locus_repair=true locus_min_support=2 --label "$T syn real constrained+veto+r3+locus(1sp,min2)" &
wait
python scripts/ob_merge_iter.py $R --arm default --synteny shuffled --syn shuffle_seed=11 $C locus_repair=true locus_min_support=2 --label "$T syn shuffled(s11) constrained+veto+r3+locus(1sp,min2)"
