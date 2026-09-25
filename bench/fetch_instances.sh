#!/usr/bin/env bash
# Downloads public benchmark instances into bench/instances/<set>/ as plain .mps files.
#
#   bench/fetch_instances.sh netlib           # Netlib LP (decoded from netlib's compressed EMPS format)
#   bench/fetch_instances.sh netlib-infeas    # Netlib infeasible LPs
#   bench/fetch_instances.sh miplib           # MIPLIB 2017 benchmark set
#
# Instances are never committed (see .gitignore). Requires curl, a C compiler, gunzip and unzip.
set -euo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
dest_root="${here}/instances"

netlib_names=(
  25fv47 80bau3b adlittle afiro agg agg2 agg3 bandm beaconfd blend bnl1 bnl2 boeing1 boeing2
  bore3d brandy capri cycle czprob d2q06c d6cube degen2 degen3 dfl001 e226 etamacro fffff800
  finnis fit1d fit1p fit2d fit2p forplan ganges gfrd-pnc greenbea greenbeb grow15 grow22 grow7
  israel kb2 lotfi maros maros-r7 modszk1 nesm perold pilot pilot.ja pilot.we pilot4 pilot87
  pilotnov recipe sc105 sc205 sc50a sc50b scagr25 scagr7 scfxm1 scfxm2 scfxm3
  scorpion scrs8 scsd1 scsd6 scsd8 sctap1 sctap2 sctap3 seba share1b share2b shell ship04l
  ship04s ship08l ship08s ship12l ship12s sierra stair standata standgub standmps stocfor1
  stocfor2 tuff vtp.base wood1p woodw
)

netlib_infeas_names=(
  bgdbg1 bgetam bgindy bgprtr box1 ceria3d chemcom cplex1 cplex2 ex72a ex73a forest6 galenet
  gosh gran greenbea itest2 itest6 klein1 klein2 klein3 mondou2 pang pilot4i qual reactor
  refinery vol1 woodinfe
)

# fetch_emps <netlib directory> <destination> <names...>
fetch_emps() {
  local dir="$1" dest="$2"
  shift 2
  mkdir -p "${dest}"
  local emps="${dest}/.emps"
  if [[ ! -x "${emps}" ]]; then
    curl -fsSL https://www.netlib.org/lp/data/emps.c -o "${dest}/.emps.c"
    cc -O2 -o "${emps}" "${dest}/.emps.c"
  fi
  for name in "$@"; do
    local out="${dest}/${name}.mps"
    [[ -s "${out}" ]] && continue
    echo "netlib: ${name}"
    if curl -fsSL "https://www.netlib.org/lp/${dir}/${name}" | "${emps}" > "${out}.tmp" \
        && [[ -s "${out}.tmp" ]]; then
      mv "${out}.tmp" "${out}"
    else
      rm -f "${out}.tmp"
      echo "netlib: ${name} not available, skipped" >&2
    fi
  done
}

fetch_miplib() {
  local dest="${dest_root}/miplib2017"
  mkdir -p "${dest}"
  local zip="${dest}/.benchmark.zip"
  [[ -s "${zip}" ]] || curl -fSL https://miplib.zib.de/downloads/benchmark.zip -o "${zip}"
  unzip -oq "${zip}" -d "${dest}"
  find "${dest}" -name '*.mps.gz' -exec gunzip -f {} +
  curl -fsSL https://miplib.zib.de/downloads/miplib2017-v31.solu -o "${dest}/miplib2017.solu" || true
}

if [[ $# -eq 0 ]]; then
  echo "usage: $0 netlib|netlib-infeas|miplib ..." >&2
  exit 2
fi
for set in "$@"; do
  case "${set}" in
    netlib) fetch_emps data "${dest_root}/netlib" "${netlib_names[@]}" ;;
    netlib-infeas) fetch_emps infeas "${dest_root}/netlib-infeas" "${netlib_infeas_names[@]}" ;;
    miplib) fetch_miplib ;;
    *) echo "unknown instance set '${set}'" >&2; exit 2 ;;
  esac
done
