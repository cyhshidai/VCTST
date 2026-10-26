#!/usr/bin/env bash
# Reproduce the four prune log blocks of README section 9 (no DFT code needed).
#   bash examples/prune_demo.sh
# The generated *.yaml are written next to examples/lj_ar.yaml (the endpoint
# paths inside it are "./structures/...", which vcneb resolves relative to the
# YAML file, not to your current directory) and are removed again at the end.
set -u
cd "$(dirname "$0")/.." || exit 1
base=examples/lj_ar.yaml
make_case () {  # $1=file  $2=min_distance  $3=prune
    python3 - "$base" "$1" "$2" "$3" <<'PY'
import sys
src, dst, dmin, prune = sys.argv[1:5]
s = open(src).read()
s = s.replace("  images: 5", "  images: 8")
s = s.replace("  max_steps: 400",
              "  max_steps: 400\n  min_distance: %s\n  prune: %s" % (dmin, prune))
s = s.replace("./runs/lj_ar", "./runs/pr_%s_%s" % (dmin, prune))
open(dst, "w").write(s)
PY
}

make_case examples/pr_31a.yaml 3.1 auto
make_case examples/pr_33a.yaml 3.3 auto
make_case examples/pr_35e.yaml 3.5 error
make_case examples/pr_35a.yaml 3.5 auto

for f in 31a 33a 35e 35a; do
    echo "================ examples/pr_$f.yaml"
    vcneb run "examples/pr_$f.yaml" --dry-run 2>&1 | grep -v '^$'
    echo
done

rm -f examples/pr_31a.yaml examples/pr_33a.yaml examples/pr_35e.yaml \
      examples/pr_35a.yaml
echo "all four cases above are README section 9; runs/ and the temp yaml are ignored by git"
