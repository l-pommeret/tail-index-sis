#!/bin/sh
# Commit and push every finished grid cell as soon as it lands on disk.
# A cell file is only staged once it has been untouched for 20 seconds, so a
# half-written .rds is never committed.  Exits when all cells are in, or when
# the driver is gone and nothing is left to pick up.
# usage: code/sh/git_watch_cells.sh [CELLDIR] [TOTAL]
set -u
DIR=${1:-results/grid/comparison_cells}
TOTAL=${2:-288}
export GIT_TERMINAL_PROMPT=0

while :; do
  tracked=$(git ls-files "$DIR")
  new=0
  for f in $(find "$DIR" -name '*.rds' ! -newermt '-20 seconds' 2>/dev/null | sort); do
    case "
$tracked" in *"
$f"*) continue ;; esac
    git add "$f" || continue
    if git commit -q -m "grid cell $(basename "$f" .rds): 40 replications"; then
      new=$((new + 1))
    fi
  done
  if [ "$new" -gt 0 ]; then
    # refresh the human-readable report from every cell available so far
    if PATH="$PWD/.renv/bin:$PATH" KAPPA=0.20 ASTAR=0.30 BSTAR=0.10 \
       GRID_NCELL="$TOTAL" Rscript code/R/summarize_grid.R "$DIR" \
       results/grid >/dev/null 2>&1; then
      git add results/grid/rapport_grille.md results/grid/grid_summary.csv \
              results/grid/grid_timings.csv
      git commit -q -m "grid: refresh report ($(git ls-files "$DIR" | wc -l)/$TOTAL cells)" \
        >/dev/null 2>&1 || true
    fi
    git push -q origin HEAD 2>&1 | tail -2
    echo "$(date +%H:%M:%S) pushed $new cell(s), $(git ls-files "$DIR" | wc -l)/$TOTAL total"
  fi
  committed=$(git ls-files "$DIR" | wc -l)
  if [ "$committed" -ge "$TOTAL" ]; then
    echo "all $TOTAL cells committed and pushed"; break
  fi
  if ! pgrep -f "[r]un_grid.R" >/dev/null 2>&1; then
    if [ "$(find "$DIR" -name '*.rds' | wc -l)" -le "$committed" ]; then
      echo "driver finished; $committed cells committed and pushed"; break
    fi
  fi
  sleep 30
done
