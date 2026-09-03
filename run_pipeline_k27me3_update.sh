#!/bin/bash
set -e

echo "Starting step 7d..."
out_7d=$(python3 Scripts/step_7d_classify_regions.py)
echo "$out_7d"
id_7d=$(echo "$out_7d" | grep -oP 'Job <\K\d+' || true)

if [ -n "$id_7d" ]; then
    echo "Waiting for step 7d (Jobs: $id_7d)..."
    bash Scripts/utils/poll_jobs.sh step_7d $id_7d --timeout 7200
else
    echo "Failed to get Job ID for step 7d!"
    exit 1
fi

echo "Starting step 8a for K27me3..."
out_8a_1=$(python3 Scripts/step_8a_computematrix.py --base_prefix retained_h3k27me3 --regions Analysis_Data/k27me3_classification/retained_k27me3.bed)
echo "$out_8a_1"
out_8a_2=$(python3 Scripts/step_8a_computematrix.py --base_prefix lost_h3k27me3 --regions Analysis_Data/k27me3_classification/lost_k27me3.bed)
echo "$out_8a_2"
ids_8a=$(echo -e "$out_8a_1\n$out_8a_2" | grep -oP 'Job <\K\d+' | tr '\n' ' ' || true)

if [ -n "$ids_8a" ]; then
    echo "Waiting for step 8a (Jobs: $ids_8a)..."
    bash Scripts/utils/poll_jobs.sh step_8a $ids_8a --timeout 7200
else
    echo "Failed to get Job IDs for step 8a!"
    exit 1
fi

echo "Starting step 8b for all 4..."
out_8b_1=$(python3 Scripts/step_8b_plotheatmap.py --base_prefix retained_h3k27me2)
echo "$out_8b_1"
out_8b_2=$(python3 Scripts/step_8b_plotheatmap.py --base_prefix lost_h3k27me2)
echo "$out_8b_2"
out_8b_3=$(python3 Scripts/step_8b_plotheatmap.py --base_prefix retained_h3k27me3)
echo "$out_8b_3"
out_8b_4=$(python3 Scripts/step_8b_plotheatmap.py --base_prefix lost_h3k27me3)
echo "$out_8b_4"
ids_8b=$(echo -e "$out_8b_1\n$out_8b_2\n$out_8b_3\n$out_8b_4" | grep -oP 'Job <\K\d+' | tr '\n' ' ' || true)

if [ -n "$ids_8b" ]; then
    echo "Waiting for step 8b (Jobs: $ids_8b)..."
    bash Scripts/utils/poll_jobs.sh step_8b $ids_8b --timeout 7200
else
    echo "Failed to get Job IDs for step 8b!"
    exit 1
fi

echo "Pipeline update finished!"
