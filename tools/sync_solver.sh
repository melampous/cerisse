#!/bin/bash
# Push the local solver source to the deployment mirrors and verify.
#
# The local git repository is the only place the solver is edited. AWS and
# Y9000X hold deploy copies of src/ alone: they are not git repositories, and
# editing them is what produced the 25 and 23 divergent copies that this
# arrangement replaces. Each machine keeps its own lib/amrex, because the two
# are pinned by different compiler and CUDA constraints (AWS needs the
# CCCL-guarded AMReX to build against CUDA 12.6).
#
# After rsync the script compares an md5 manifest of every file under src/, so
# a partial transfer or a stale file is reported rather than assumed away.
#
#   ./tools/sync_solver.sh            # both mirrors
#   ./tools/sync_solver.sh aws        # one of: aws, y9000x
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
AWS_DEST=/shared/cerisse_unified_20260808
Y_DEST=/home/melampous/cerisse_unified_20260808
TARGETS="${1:-all}"

cd "$REPO" || exit 1

dirty=$(git status --porcelain -- src/ | wc -l)
branch=$(git rev-parse --abbrev-ref HEAD)
head=$(git rev-parse --short HEAD)
echo "local: $branch @ $head   ($dirty uncommitted change(s) under src/)"
if [ "$dirty" -ne 0 ]; then
  echo "  note: uncommitted work will be pushed to the mirrors but is not in git."
fi

manifest=$(mktemp)
(cd src && find . -type f | sort | xargs md5sum | sed 's|  ./|  |') > "$manifest"
echo "manifest: $(wc -l < "$manifest") files"

verify() {  # verify <ssh-target> <remote-src-dir> <label>
  scp -q "$manifest" "$1:/tmp/solver_manifest.txt" || { echo "  $3: manifest copy FAILED"; return 1; }
  ssh "$1" "cd $2 && find . -type f | sort | xargs md5sum | sed 's|  ./|  |' > /tmp/solver_remote.txt
    d=\$(diff <(sort /tmp/solver_manifest.txt) <(sort /tmp/solver_remote.txt) | wc -l)
    if [ \"\$d\" -eq 0 ]; then echo '  $3: OK, byte-for-byte identical'
    else echo \"  $3: MISMATCH (\$d lines)\"; diff <(sort /tmp/solver_manifest.txt) <(sort /tmp/solver_remote.txt) | head -10; fi"
}

if [ "$TARGETS" = all ] || [ "$TARGETS" = aws ]; then
  echo "-> AWS $AWS_DEST"
  rsync -a --delete src/ "aws-hpc:$AWS_DEST/src/" && verify aws-hpc "$AWS_DEST/src" AWS
fi

if [ "$TARGETS" = all ] || [ "$TARGETS" = y9000x ]; then
  echo "-> Y9000X $Y_DEST"
  rsync -a --delete --chmod=u+w src/ "Y9000X:$Y_DEST/src/" && verify Y9000X "$Y_DEST/src" Y9000X
fi

rm -f "$manifest"
