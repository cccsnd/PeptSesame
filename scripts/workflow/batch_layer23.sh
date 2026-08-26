#!/usr/bin/env bash
# Layer2+3 batch execution (serial, per-species L2->L3->self-check)
# Usage: bash batch_layer23.sh <species1> <species2> ...
set -u
cd "$(dirname "$0")/../.."   # repo root
PY=${PYTHON:-python}
R2=results

for sp in "$@"; do
    echo "[$(date +%H:%M:%S)] === $sp Layer2 ==="
    $PY scripts/workflow/run_layer2.py "$sp" > $R2/02_layer2_scoring/${sp}_layer2.log 2>&1
    L2_OK=$(tail -1 $R2/02_layer2_scoring/${sp}_layer2.log | grep -c "自检通过")
    if [ "$L2_OK" != "1" ]; then
        echo "[$sp] Layer2 失败! 日志尾部:"; tail -5 $R2/02_layer2_scoring/${sp}_layer2.log; continue
    fi
    echo "[$(date +%H:%M:%S)] === $sp Layer3 ==="
    $PY scripts/workflow/run_layer3.py "$sp" > $R2/03_layer3_classify/${sp}_layer3.log 2>&1
    L3_OK=$(tail -1 $R2/03_layer3_classify/${sp}_layer3.log | grep -c "自检通过")
    if [ "$L3_OK" != "1" ]; then
        echo "[$sp] Layer3 失败! 日志尾部:"; tail -5 $R2/03_layer3_classify/${sp}_layer3.log; continue
    fi
    echo "[$(date +%H:%M:%S)] $sp 完成 ✅ (L2+L3)"
done
echo "[$(date +%H:%M:%S)] 全部批次结束"
