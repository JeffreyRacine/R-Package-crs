#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
ATTR_DIR="${ROOT_DIR}/src/nomad4_src/src/Attribute"
OUT_DIR="${ROOT_DIR}/inst/nomad"
OUT_FILE="${OUT_DIR}/NOMAD_4_5_0_OPTIONS_REFERENCE.md"

mkdir -p "${OUT_DIR}"

declare -a FILE_ORDER=(
  "pbAttributesDefinition.txt"
  "runAttributesDefinition.txt"
  "runAttributesDefinitionCOOP.txt"
  "runAttributesDefinitionCS.txt"
  "runAttributesDefinitionDMulti.txt"
  "runAttributesDefinitionDisco.txt"
  "runAttributesDefinitionIBEX.txt"
  "runAttributesDefinitionLH.txt"
  "runAttributesDefinitionNM.txt"
  "runAttributesDefinitionPSDSSD.txt"
  "runAttributesDefinitionQPSolver.txt"
  "runAttributesDefinitionQuadModel.txt"
  "runAttributesDefinitionSgtelibModel.txt"
  "runAttributesDefinitionVNS.txt"
  "evalAttributesDefinition.txt"
  "displayAttributesDefinition.txt"
  "cacheAttributesDefinition.txt"
  "evaluatorControlAttributesDefinition.txt"
  "evaluatorControlGlobalAttributesDefinition.txt"
)

title_for_file() {
  case "$1" in
    pbAttributesDefinition.txt) echo "Problem Parameters" ;;
    runAttributesDefinition.txt) echo "Run Parameters (Core)" ;;
    runAttributesDefinitionCOOP.txt) echo "Run Parameters (COOP)" ;;
    runAttributesDefinitionCS.txt) echo "Run Parameters (Coordinate Search)" ;;
    runAttributesDefinitionDMulti.txt) echo "Run Parameters (DMultiMads)" ;;
    runAttributesDefinitionDisco.txt) echo "Run Parameters (DiscoMads)" ;;
    runAttributesDefinitionIBEX.txt) echo "Run Parameters (IBEX)" ;;
    runAttributesDefinitionLH.txt) echo "Run Parameters (Latin Hypercube)" ;;
    runAttributesDefinitionNM.txt) echo "Run Parameters (Nelder-Mead)" ;;
    runAttributesDefinitionPSDSSD.txt) echo "Run Parameters (PSD-Mads)" ;;
    runAttributesDefinitionQPSolver.txt) echo "Run Parameters (QP Solver)" ;;
    runAttributesDefinitionQuadModel.txt) echo "Run Parameters (Quadratic Model)" ;;
    runAttributesDefinitionSgtelibModel.txt) echo "Run Parameters (Sgtelib Model)" ;;
    runAttributesDefinitionVNS.txt) echo "Run Parameters (VNS)" ;;
    evalAttributesDefinition.txt) echo "Evaluation Parameters" ;;
    displayAttributesDefinition.txt) echo "Display Parameters" ;;
    cacheAttributesDefinition.txt) echo "Cache Parameters" ;;
    evaluatorControlAttributesDefinition.txt) echo "Evaluator Control Parameters" ;;
    evaluatorControlGlobalAttributesDefinition.txt) echo "Evaluator Control Global Parameters" ;;
    *) echo "$1" ;;
  esac
}

extract_rows() {
  awk '
function trim(s) {
  sub(/^[[:space:]]+/, "", s)
  sub(/[[:space:]]+$/, "", s)
  return s
}
function esc_md(s) {
  gsub(/\|/, "\\|", s)
  return s
}
{
  if ($0 ~ /^[A-Za-z][A-Za-z0-9_]*$/) {
    name = $0
    if ((getline type) <= 0) next
    if ((getline defv) <= 0) next
    short = ""
    while ((getline) > 0) {
      if ($0 ~ /^\\\(/) {
        short = $0
        break
      }
    }
    sub(/^\\\([[:space:]]*/, "", short)
    sub(/[[:space:]]*\\\)$/, "", short)
    printf("| `%s` | `%s` | `%s` | %s |\n",
           esc_md(name), esc_md(trim(type)), esc_md(trim(defv)), esc_md(trim(short)))
  }
}
'
}

{
  printf "%s\n\n" "# NOMAD 4.5.0 Option Reference for \`crs\`"
  printf "Generated on: %s\n\n" "$(date '+%Y-%m-%d')"
  printf "Source of truth: vendored NOMAD definition files in \`%s\`\n\n" "${ATTR_DIR#"${ROOT_DIR}"/}"
  printf "%s\n" "Scope:"
  printf "%s\n" "- This reference is generated from embedded NOMAD 4.5.0 sources used by \`crs\`."
  printf "%s\n" "- \`snomadr(opts=...)\` passes option names/values through to NOMAD."
  printf "%s\n" "- \`display.nomad.progress=FALSE\` also sets \`DISPLAY_DEGREE 0\`, \`DISPLAY_ALL_EVAL false\`, and \`DISPLAY_UNSUCCESSFUL false\`."
  printf "%s\n" "- \`nomad.opt\` in the working directory is read after \`opts\` and can override earlier values."
  printf "%s\n\n" "- Option names use NOMAD 4.5.0 terminology directly (for example, \`MIN_FRAME_SIZE\`)."

  total=0
  for filename in "${FILE_ORDER[@]}"; do
    path="${ATTR_DIR}/${filename}"
    if [[ ! -f "${path}" ]]; then
      continue
    fi
    rows_file="$(mktemp)"
    extract_rows < "${path}" > "${rows_file}"
    count="$(wc -l < "${rows_file}" | tr -d " ")"
    total=$((total + count))

    printf "## %s (%s)\n\n" "$(title_for_file "${filename}")" "${filename}"
    printf "Total options in this section: %s\n\n" "${count}"
    printf "| Option | Type | Default | Short Description |\n"
    printf "| --- | --- | --- | --- |\n"
    cat "${rows_file}"
    printf "\n"
    rm -f "${rows_file}"
  done

  printf "## Totals\n\n"
  printf "%s\n" "- Total option entries across all definition files: ${total}"
  printf "%s\n" "- To inspect option help at runtime: \`snomadr(information = list(\"help\" = \"-h\"))\`"
  printf "%s\n" "- To inspect a specific option at runtime: \`snomadr(information = list(\"help\" = \"-h OPTION_NAME\"))\`"
} > "${OUT_FILE}"

echo "Wrote ${OUT_FILE}"
