#! /bin/bash
set -e

# Cleanup function to handle script termination
# This function will be called on script exit or interruption
cleanup() {
    pkill -P $$  # Kill all the child processes of the current process group
    # Optional: Delete temporary files
    [ -f "$TMP_FILE" ] && rm "$TMP_FILE"
    exit 1
}

# Register Signal Capture
trap 'cleanup' INT TERM EXIT

printf "[ProType] [Size] [Metric] [Dim] [Delta] [Online_Com.(MB)] [Time(s)] [Offline_Com.(MB)] [Offline_Time(s)]\n"


for n in "${ns[@]}"; do
    for dim in "${dims[@]}"; do
        for delta in "${deltas[@]}"; do
        ./build/fpsi  -nn $n -dim $dim -delta $delta -trait 10
        done
        echo
    done
done
