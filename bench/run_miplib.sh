#!/usr/bin/env bash
# Runs an instance list in parallel shards (single-threaded by default) and checks the results
# against the MIPLIB known solutions.
#
#   bench/run_miplib.sh [SHARDS] [SECONDS] [LIST] [SAMAYA] [OUT] [THREADS]
#
#   SHARDS   parallel runs (default 4; use the number of performance cores, fewer if RAM < 4 GB
#            per shard)
#   SECONDS  time limit per instance (default 600; 60 for quick screening)
#   LIST     instance list (default bench/miplib_small.test, fetched with
#            bench/fetch_instances.sh miplib-list bench/miplib_small.test)
#   SAMAYA   solver binary (default build/release/apps/cli/samaya)
#   OUT      output directory for the per-shard CSV and logs (default bench/results/miplib)
#   THREADS  tree-search threads per run (default 1); keep SHARDS x THREADS <= cores
#
# Prints the per-shard summaries; bench/report.py OUT/shard*.csv makes a table.
set -euo pipefail

cd "$(dirname "$0")/.."
shards="${1:-4}"
seconds="${2:-600}"
list="${3:-bench/miplib_small.test}"
samaya="${4:-build/release/apps/cli/samaya}"
out="${5:-bench/results/miplib}"
threads="${6:-1}"
dir="bench/instances/$(basename "${list%.*}")"

[[ -x "${samaya}" ]] || { echo "solver not found: ${samaya}" >&2; exit 1; }
[[ -d "${dir}" ]] || { echo "instances not found in ${dir}; run bench/fetch_instances.sh miplib-list ${list}" >&2; exit 1; }
mkdir -p "${out}"

mapfile -t files < <(grep -v '^#' "${list}" | grep -v '^$' | sed "s|^|${dir}/|; s|$|.mps|")
for ((k = 0; k < shards; ++k)); do
  shard=()
  for ((i = k; i < ${#files[@]}; i += shards)); do shard+=("${files[i]}"); done
  [[ ${#shard[@]} -gt 0 ]] || continue
  python3 bench/harness.py "${shard[@]}" --samaya "${samaya}" --time-limit "${seconds}" \
    --threads "${threads}" --solu "${dir}/miplib2017.solu" --out "${out}/shard${k}.csv" > "${out}/shard${k}.log" 2>&1 &
done
wait
cat "${out}"/shard*.log | grep -E "Known solutions|WRONG|solved"
python3 - "${out}" <<'EOF'
import csv, glob, sys
rows = [r for f in glob.glob(sys.argv[1] + "/shard*.csv") for r in csv.DictReader(open(f))]
solved = sum(r["status"] in ("optimal", "infeasible") for r in rows)
print(f"total: {solved}/{len(rows)} solved (optimal or proven infeasible)")
EOF
