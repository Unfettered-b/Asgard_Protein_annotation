#!/bin/bash

# Script to repeatedly run esm_inference.py until it succeeds

while true; do
    echo "Starting esm_inference.py..."
    
    python scripts/esm_inference.py
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo "esm_inference.py completed successfully!"
        break
    else
        echo "esm_inference.py failed with exit code $exit_code. Retrying..."
        sleep 2   # small delay before retry (optional)
    fi
done
