#!/bin/bash

# Define the input and output file paths
input_file="/Users/leobailloeul/Documents/coding/decaysimulation/mesongen/output-data/mesons_seed1234bigjsi.root"
output_file="/Users/leobailloeul/Documents/coding/decaysimulation/decay/output-data/mcp-production-decay-vector-meson.root"
efficiency_output="/Users/leobailloeul/Documents/coding/decaysimulation/decay/sensitivity-plot/efficiency_output.txt"
total_efficiency_output="/Users/leobailloeul/Documents/coding/decaysimulation/decay/sensitivity-plot/total_efficiency_output.txt"

# Ensure the total efficiency output file is empty before appending new data
> $total_efficiency_output

# Read mchi values from mchi_values.txt
while IFS= read -r mchi; do
    echo "Running decayVectorMeson.cc with mchi = $mchi"

    # Clear the efficiency output file before each run
    > $efficiency_output

    # Run the C++ script with the current mchi value using the relative path
    ../decayVectorMeson "$input_file" "$output_file" "$mchi"

    # Check if the efficiency output file exists and is not empty
    if [ -s $efficiency_output ]; then
        echo "Appending efficiency to $total_efficiency_output"
        # Append the efficiency to the total output file
        cat $efficiency_output >> $total_efficiency_output
    else
        echo "Warning: $efficiency_output is empty or does not exist. Skipping."
    fi

done < mchi_values.txt

echo "All simulations completed. Total efficiencies are stored in $total_efficiency_output"