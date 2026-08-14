#!/bin/sh
# Commit every finished grid cell as soon as it lands on disk, and push in
# groups of PUSH_EVERY cells (the Markdown report is regenerated at each push,
# not at each cell, to keep the refresh off the critical path).
# A cell file is only staged once it has been untouched for 20 seconds, so a
# half-written .rds is never committed.  Exits when all cells are in, or when
# the driver is gone and nothing is left to pick up; a partial group is always
# pushed before exiting.
# usage: code/sh/git_watch_cells.sh [CELLDIR] [TOTAL] [PUSH_EVERY]
set -u
DIR=${1:-results/grid/comparison_cells}
TOTAL=${2:-288}
PUSH_EVERY=${3:-5}
export GIT_TERMINAL_PROMPT=0

publish() {   # regenerate the report, commit it, push everything pending
  committed=$(git ls-files "$DIR" | wc -l)
  if PATH="$PWD/.renv/bin:$PATH" KAPPA=0.20 ASTAR=0.30 BSTAR=0.10 \
     GRID_NCELL="$TOTAL" Rscript code/R/summarize_grid.R "$DIR" \
     results/grid >/dev/null 2>&1; then
    git add results/grid/rapport_grille.md results/grid/grid_summary.csv \
            results/grid/grid_timings.csv
    git commit -q -m "grid: refresh report ($committed/$TOTAL cells)" \
      >/dev/null 2>&1 || true
  fi
  git push -q origin HEAD 2>&1 | tail -2
  echo "$(date +%H:%M:%S) pushed, $committed/$TOTAL cells"
}

pending=0
while :; do
  tracked=$(git ls-files "$DIR")
  for f in $(find "$DIR" -name '*.rds' ! -newermt '-20 seconds' 2>/dev/null | sort); do
    case "
$tracked" in *"
$f"*) continue ;; esac
    git add "$f" || continue
    if git commit -q -m "grid cell $(basename "$f" .rds): 40 replications"; then
      pending=$((pending + 1))
    fi
  done
  if [ "$pending" -ge "$PUSH_EVERY" ]; then
    publish
    pending=0
  fi
  if [ "$(git ls-files "$DIR" | wc -l)" -ge "$TOTAL" ]; then
    [ "$pending" -gt 0 ] && publish
    echo "all $TOTAL cells committed and pushed"; break
  fi
  if ! pgrep -f "[r]un_grid.R" >/dev/null 2>&1; then
    if [ "$(find "$DIR" -name '*.rds' | wc -l)" -le "$(git ls-files "$DIR" | wc -l)" ]; then
      [ "$pending" -gt 0 ] && publish
      echo "driver finished; $(git ls-files "$DIR" | wc -l) cells committed and pushed"
      break
    fi
  fi
  sleep 30
done
