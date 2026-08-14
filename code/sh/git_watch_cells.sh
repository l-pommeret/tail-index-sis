#!/bin/sh
# Commit every finished grid cell as soon as it lands on disk, and push once per
# driver batch: the driver writes all the cells of a batch back to back and then
# spends minutes on the next one, so a burst that has gone quiet is a batch
# boundary.  The Markdown report is regenerated at each push.
# A cell file is only staged once it has been untouched for QUIET seconds, so a
# half-written .rds is never committed, and a push only fires once nothing is
# being written any more.  A partial group is always pushed before exiting.
# usage: code/sh/git_watch_cells.sh [CELLDIR] [TOTAL] [MAXHOLD]
set -u
DIR=${1:-results/grid/comparison_cells}
TOTAL=${2:-288}
MAXHOLD=${3:-40}          # safety valve: never hold more than this many cells
QUIET=20                  # seconds of inactivity that mark a finished burst
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
  echo "$(date +%H:%M:%S) batch pushed ($pending cells), $committed/$TOTAL total"
}

pending=0
while :; do
  tracked=$(git ls-files "$DIR")
  for f in $(find "$DIR" -name '*.rds' ! -newermt "-$QUIET seconds" 2>/dev/null | sort); do
    case "
$tracked" in *"
$f"*) continue ;; esac
    git add "$f" || continue
    if git commit -q -m "grid cell $(basename "$f" .rds): 40 replications"; then
      pending=$((pending + 1))
    fi
  done
  # a batch is complete when nothing has been written for QUIET seconds
  writing=$(find "$DIR" -name '*.rds' -newermt "-$QUIET seconds" 2>/dev/null | wc -l)
  if [ "$pending" -gt 0 ] && { [ "$writing" -eq 0 ] || [ "$pending" -ge "$MAXHOLD" ]; }; then
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
