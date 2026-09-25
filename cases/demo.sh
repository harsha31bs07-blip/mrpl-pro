#!/usr/bin/env bash
# End-to-end MRPL demo: generate the case-study models, solve each with samaya (every result is
# checked by the independent verifier), and print a readable plan.
#
#   cases/demo.sh [SIZE] [SAMAYA]
#
#   SIZE     small (default), medium or large
#   SAMAYA   solver binary (default build/release/apps/cli/samaya)
set -euo pipefail

cd "$(dirname "$0")/.."
size="${1:-small}"
samaya="${2:-build/release/apps/cli/samaya}"
out="cases/instances"
if [[ ! -x "${samaya}" ]]; then
  echo "solver not found: ${samaya}" >&2
  echo "build it with: cmake --preset release && cmake --build --preset release" >&2
  exit 1
fi

python3 cases/mrpl.py generate --out "${out}" > /dev/null
for family in plan crude utility; do
  model="${out}/mrpl_${family}_${size}.mps"
  echo "=================================================================================="
  echo "${model}"
  echo "----------------------------------------------------------------------------------"
  start=$(date +%s.%N)
  "${samaya}" --time-limit 60 --log-level 0 --solution "${out}/${family}_${size}.sol" "${model}"
  end=$(date +%s.%N)
  python3 cases/mrpl.py report "${model}" "${out}/${family}_${size}.sol"
  printf '\nsolved in %.2f s\n\n' "$(echo "${end} - ${start}" | bc)"
done
