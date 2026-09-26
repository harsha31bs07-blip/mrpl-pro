#!/usr/bin/env bash
# A/B comparison of two samaya builds on the MIPLIB screening set, both builds running at the same
# time so machine load affects them equally (two instances at a time: four processes).
#
#   bench/ab_bound.sh OUT BIN_A BIN_B [SECONDS] [LIST]
#
# Writes OUT/<instance>.{a,b}.json (the CLI's --json line) and OUT/DONE; then
#   python3 bench/ab_report.py OUT
# compares the proven bounds with the known optima. Use the bound, not the incumbent, to judge a
# change to cuts or the search: the incumbent's gap moves by several points between two runs of
# the same binary.
set -euo pipefail
cd "$(dirname "$0")/.."
out="${1:?usage: bench/ab_bound.sh OUT BIN_A BIN_B [SECONDS] [LIST]}"
a="$(realpath "${2:?}")"
b="$(realpath "${3:?}")"
seconds="${4:-60}"
list="${5:-bench/miplib_small.test}"
dir="bench/instances/$(basename "${list%.*}")"
mkdir -p "${out}"
mapfile -t names < <(grep -v '^#' "${list}" | grep -v '^$')
run() {
  "$1" --threads 1 --json --log-level 0 --time-limit "${seconds}" "${dir}/$2.mps" 2>/dev/null |
    tail -1 > "${out}/$2.$3.json"
}
for ((i = 0; i < ${#names[@]}; i += 2)); do
  run "${a}" "${names[i]}" a &
  run "${b}" "${names[i]}" b &
  if ((i + 1 < ${#names[@]})); then
    run "${a}" "${names[i + 1]}" a &
    run "${b}" "${names[i + 1]}" b &
  fi
  wait
done
echo finished > "${out}/DONE"
