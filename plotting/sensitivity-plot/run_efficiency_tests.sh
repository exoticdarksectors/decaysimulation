#!/bin/bash

# Define the input and output file paths
base_directory="/Users/leobailloeul/Documents/coding/decaysimulation"

# declare -a two_body_input_files
# declare -a two_body_output_files

# two_body_input_files=(
#     #"${base_directory}/mesongen/output-data/mesons_seed1234_443.root"
#     #"${base_directory}/mesongen/output-data/mesons_seed1234_113.root"
#     #"${base_directory}/mesongen/output-data/mesons_seed1234_223.root"
#     #"${base_directory}/mesongen/output-data/mesons_seed1234_333.root"
# )
# two_body_output_files=(
#     #"${base_directory}/decay/output-data/1040mcp-production-decay-vector-meson-jsi.root"
#     #"${base_directory}/decay/output-data/1040mcp-production-decay-vector-meson-rho.root"
#     #"${base_directory}/decay/output-data/1040mcp-production-decay-vector-meson-omega.root"
#     #"${base_directory}/decay/output-data/1040mcp-production-decay-vector-meson-phi.root"
# )

# declare -a three_body_input_file
# declare -a three_body_output_file

# three_body_input_file=(
# # "${base_directory}/mesongen/output-data/mesons_seed1234_111.root"
# "${base_directory}/mesongen/output-data/mesons_seed1234_221DBG.root"
# )
# three_body_output_file=(
# # "${base_directory}/decay/output-data/mcp-production574mpi0.root"
# "${base_directory}/decay/output-data/mcp-production1040meta.root"
# )

# declare -a efficiency_output_2bodydecay
# declare -a total_efficiency_output_2bodydecay

# efficiency_output_2bodydecay=(
#     #"TwoByTwoefficiency-output2bodydecay-jsi.txt"
#     #"TwoByTwoefficiency-output2bodydecay-rho.txt"
#     #"TwoByTwoefficiency-output2bodydecay-omega.txt"
#     #"TwoByTwoefficiency-output2bodydecay-phi.txt"
# )

# total_efficiency_output_2bodydecay=(
#     #"${base_directory}/plotting/sensitivity-plot/TwoByTwototal_efficiency_output2bodydecay-jsi.txt"
#     #"${base_directory}/plotting/sensitivity-plot/TwoByTwototal_efficiency_output2bodydecay-rho.txt"
#     #"${base_directory}/plotting/sensitivity-plot/TwoByTwototal_efficiency_output2bodydecay-omega.txt"
#     #"${base_directory}/plotting/sensitivity-plot/TwoByTwototal_efficiency_output2bodydecay-phi.txt"
# )

# declare -a efficiency_output_three_body
# declare -a total_efficiency_output_three_body

# efficiency_output_three_body=(
#   # "TwoByTwoefficiency_outputdecayPion.txt"
#   "TwoByTwoefficiency_outputdecayEta.txt"
# )

# total_efficiency_output_three_body=(
#   # "${base_directory}/plotting/sensitivity-plot/TwoByTwototal_efficiency_outputdecayPion.txt"
#   "${base_directory}/plotting/sensitivity-plot/TwoByTwototal_efficiency_outputdecayEtaDBG.txt"
# )
# declare -a three-body-masses
# three_body_masses=(
# # "0.1349768"
# "0.547862"
# )


declare -a two_body_input_files
declare -a two_body_output_files

two_body_input_files=(
   "${base_directory}/mesongen/output-data/mesons_seed1234_443.root"
#    "${base_directory}/mesongen/output-data/mesons_seed1234_113.root"
#    "${base_directory}/mesongen/output-data/mesons_seed1234_223.root"
#    "${base_directory}/mesongen/output-data/mesons_seed1234_333.root"
)
two_body_output_files=(
   "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-jsi.root"
   # "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-rho.root"
   # "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-omega.root"
   # "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-phi.root"
)

declare -a three_body_input_file
declare -a three_body_output_file

three_body_input_file=(
"${base_directory}/mesongen/output-data/mesons_seed1234_111.root"
"${base_directory}/mesongen/output-data/mesons_seed1234_221.root"
)
three_body_output_file=(
"${base_directory}/decay/output-data/mcp-production574mjsi.root"
"${base_directory}/decay/output-data/mcp-production574meta.root"
)

declare -a efficiency_output_2bodydecay
declare -a total_efficiency_output_2bodydecay

efficiency_output_2bodydecay=(
   "efficiency-output2bodydecay-jsi.txt"
   # "efficiency-output2bodydecay-rho.txt"
   # "efficiency-output2bodydecay-omega.txt"
   # "efficiency-output2bodydecay-phi.txt"
)

total_efficiency_output_2bodydecay=(
   "${base_directory}/plotting/sensitivity-plot/total_efficiency_output2bodydecay-jsi574mSpin.txt"
   # "${base_directory}/plotting/sensitivity-plot/total_efficiency_output2bodydecay-rho574mSpin.txt"
   # "${base_directory}/plotting/sensitivity-plot/total_efficiency_output2bodydecay-omega574mSpin.txt"
   # "${base_directory}/plotting/sensitivity-plot/total_efficiency_output2bodydecay-phi574mSpin.txt"
)

declare -a efficiency_output_three_body
declare -a total_efficiency_output_three_body

efficiency_output_three_body=(
 "efficiency_outputdecayPion.txt"
 "efficiency_outputdecayEta.txt"
)

total_efficiency_output_three_body=(
 "${base_directory}/plotting/sensitivity-plot/total_efficiency_outputdecayPion574mSpin.txt"
 "${base_directory}/plotting/sensitivity-plot/total_efficiency_outputdecayEta574mSpin.txt"
)
declare -a three_body_masses
three_body_masses=(
"0.1349768"
"0.547862"
)


# Iterate over the array indices
for i in "${!two_body_input_files[@]}"; do

    # Ensure the total efficiency output file is empty before appending new data
    > "${total_efficiency_output_2bodydecay[$i]}"

    # Read mchi values from mchi_values.txt
    while IFS= read -r mchi; do
        echo "Running decayVectorMeson with mchi = $mchi"

        # Clear the efficiency output file before each run
        > "${efficiency_output_2bodydecay[$i]}"

        # Run the C++ script with the current mchi value using the relative path
        command="/Users/leobailloeul/Documents/coding/decaysimulation/decay/cmake-build-release/decayVectorMeson ${two_body_input_files[$i]} ${two_body_output_files[$i]} $mchi ${efficiency_output_2bodydecay[$i]}"
        echo "Executing: $command"
        $command

        # Check the return status of the command
        if [ $? -ne 0 ]; then
            echo "Error: decayVectorMeson failed for mchi = $mchi"
            continue
        fi

        # Check if the efficiency output file exists and is not empty
        if [ -s "${efficiency_output_2bodydecay[$i]}" ]; then
            echo "Appending efficiency to ${total_efficiency_output_2bodydecay[$i]}"
            # Append the efficiency to the total output file
            cat "${efficiency_output_2bodydecay[$i]}" >> "${total_efficiency_output_2bodydecay[$i]}"
        else
            echo "Warning: ${efficiency_output_2bodydecay[$i]} is empty or does not exist. Skipping."
        fi

    done < mchi_values.txt

done

echo "All simulations completed. Total efficiencies are stored in ${total_efficiency_output_2bodydecay[@]}"

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
    command="/Users/leobailloeul/Documents/coding/decaysimulation/decay/cmake-build-release/decayPion ${three_body_input_file[$k]} ${three_body_output_file[$k]} $mchi ${three_body_masses[$k]} ${efficiency_output_three_body[$k]}"
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