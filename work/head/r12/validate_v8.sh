#!/bin/bash
# Validate referee v8 (staging copy): head solution with matchbal, untouched round-11 code, matchbal without chunks.
cd /d/ETH-loop/work/head/r12 || exit 1
PY=/d/ETH-loop/.venv/Scripts/python.exe
run() {
  mkdir -p val/$1
  $PY -I referee_v8/referee.py "$2" val/$1/referee.json > val/$1/referee.out 2>&1
  echo "$(date +%H:%M:%S) done $1" >> val/progress.log
}
mkdir -p val
echo "$(date +%H:%M:%S) start" > val/progress.log
run good /d/ETH-loop/work/head/r12/hs12 &
run current /d/ETH-loop/algo_appd &
run bad /d/ETH-loop/work/head/r12/hs12_bad &
wait
echo "$(date +%H:%M:%S) all done" >> val/progress.log
