#!/bin/bash

# Define the input and output file paths
base_directory="/Users/leobailloeul/Documents/coding/decaysimulation"

declare -a three_body_input_file
declare -a three_body_output_file

three_body_input_file=(
"${base_directory}/mesongen/output-data/mesons_seed1234_eta.root"
)
three_body_output_file=(
"${base_directory}/decay/output-data/mcp-production40meta.root"
)

declare -a efficiency_output_three_body
declare -a total_efficiency_output_three_body

efficiency_output_three_body=(
  "efficiency_outputdecayEta.txt"
)

total_efficiency_output_three_body=(
  "${base_directory}/decay/sensitivity-plot/total_efficiency_outputdecayEta.txt"
)

declare -a three-body-masses
three_body_masses=(
"0.547862"
)

# Iterate over the array indices
for k in "${!three_body_input_file[@]}"; do

# Ensure the total efficiency output file is empty before appending new data
> "${total_efficiency_output_three_body[$k]}"

# Read mchi values from mchi_values.txt
while IFS= read -r mchi; do
    echo "Running decayPion with mchi = $mchi"

    # Clear the efficiency output file before each run
    > "${efficiency_output_three_body[$k]}"

    # Run the C++ script with the current mchi value using the relative path
    command="../cmake-build-debug/decayPion ${three_body_input_file[$k]} ${three_body_output_file[$k]} $mchi ${three_body_masses[$k]} ${efficiency_output_three_body[$k]}"
            echo "Executing: $command"
            $command

    # Check the return status of the command
    if [ $? -ne 0 ]; then
        echo "Error: decayPion failed for mchi = $mchi"
        continue
    fi

   # Check if the efficiency output file exists and is not empty
          if [ -s "${efficiency_output_three_body[$k]}" ]; then
              echo "Appending efficiency to ${total_efficiency_output_three_body[$k]}"
              # Append the efficiency to the total output file
              cat "${efficiency_output_three_body[$k]}" >> "${total_efficiency_output_three_body[$k]}"
          else
              echo "Warning: ${efficiency_output_three_body[$k]} is empty or does not exist. Skipping."
          fi

      done < mchi_values.txt

  done

echo "All simulations completed. Total efficiencies are stored in ${total_efficiency_output_three_body[@]}"