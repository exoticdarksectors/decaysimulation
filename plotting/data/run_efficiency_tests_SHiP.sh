#!/bin/bash

# Define the input and output file paths
base_directory="/Users/leobailloeul/Documents/coding/decaysimulation"
base_directory_archive="/Users/leobailloeul/Documents/coding/Archive-dec-9-2025"

declare -a two_body_input_files
declare -a two_body_output_files

two_body_input_files=(
   "${base_directory}/mesongen/output-data/mesons_seed1234_443_SHIP.root"
   "${base_directory}/mesongen/output-data/mesons_seed1234_113_SHIP.root"
   "${base_directory}/mesongen/output-data/mesons_seed1234_223_SHIP.root"
   "${base_directory}/mesongen/output-data/mesons_seed1234_333_SHIP.root"
)
two_body_output_files=(
   "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-jsi-SHIP-4by4.root"
   "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-rho-SHIP-4by4.root"
   "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-omega-SHIP-4by4.root"
   "${base_directory}/decay/output-data/mcp-production-decay-vector-meson-phi-SHIP-4by4.root"
)

declare -a three_body_input_file
declare -a three_body_output_file

three_body_input_file=(
"${base_directory}/mesongen/output-data/mesons_seed1234_111_SHIP.root"
"${base_directory}/mesongen/output-data/mesons_seed1234_221_SHIP.root"
)
three_body_output_file=(
"${base_directory}/decay/output-data/mcp-production-pi0-4by4.root"
"${base_directory}/decay/output-data/mcp-production-eta-4by4.root"
)

declare -a efficiency_output_2bodydecay
declare -a total_efficiency_output_2bodydecay

efficiency_output_2bodydecay=(
   "efficiency-output2bodydecay-jsi-4by4.txt"
   "efficiency-output2bodydecay-rho-4by4.txt"
   "efficiency-output2bodydecay-omega-4by4.txt"
   "efficiency-output2bodydecay-phi-4by4.txt"

)

total_efficiency_output_2bodydecay=(
    # "${base_directory_archive}/data-backup/total_efficiency_output2body_decay-jsi-SHiP4by4.txt"
    # "${base_directory_archive}/data-backup/total_efficiency_output2bodydecay-rho-SHiP4by4.txt"
    # "${base_directory_archive}/data-backup/total_efficiency_output2bodydecay-omega-SHiP4by4.txt"
    # "${base_directory_archive}/data-backup/total_efficiency_output2bodydecay-phi-SHiP4by4.txt"
    "${base_directory_archive}/data-backup/total_efficiency_output2body_decay-jsi-SHiP.txt"
    "${base_directory_archive}/data-backup/total_efficiency_output2bodydecay-rho-SHiP.txt"
    "${base_directory_archive}/data-backup/total_efficiency_output2bodydecay-omega-SHiP.txt"
    "${base_directory_archive}/data-backup/total_efficiency_output2bodydecay-phi-SHiP.txt"
)

declare -a efficiency_output_three_body
declare -a total_efficiency_output_three_body

efficiency_output_three_body=(
     "efficiency_outputdecayPion-4by4.txt"
     "efficiency_outputdecayPion-4by4.txt"
)

total_efficiency_output_three_body=(
  "${base_directory_archive}/data-backup/total_efficiency_outputdecayPion-SHiP.txt"
  "${base_directory_archive}/data-backup/total_efficiency_outputdecayEta-SHiP.txt"
  # "${base_directory_archive}/data-backup/total_efficiency_outputdecayPion-SHiP4by4.txt"
  # "${base_directory_archive}/data-backup/total_efficiency_outputdecayEta-SHiP4by4.txt"
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
        command="/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/decay-backup/build/decayVectorMeson ${two_body_input_files[$i]} ${two_body_output_files[$i]} $mchi ${efficiency_output_2bodydecay[$i]}"
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

    done < mship_values.txt

done

echo "All simulations completed. Total efficiencies are stored in ${total_efficiency_output_2bodydecay[@]}"

# # Iterate over the array indices
# for k in "${!three_body_input_file[@]}"; do

# # Ensure the total efficiency output file is empty before appending new data
# > "${total_efficiency_output_three_body[$k]}"

# # Read mchi values from mchi_values.txt
# while IFS= read -r mchi; do
#     echo "Running decayPion with mchi = $mchi"

#     # Clear the efficiency output file before each run
#     > "${efficiency_output_three_body[$k]}"

#     # Run the C++ script with the current mchi value using the relative path
#     command="/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/decay-backup/build/decayPion ${three_body_input_file[$k]} ${three_body_output_file[$k]} $mchi ${three_body_masses[$k]} ${efficiency_output_three_body[$k]}"
#             echo "Executing: $command"
#             $command

#     # Check the return status of the command
#     if [ $? -ne 0 ]; then
#         echo "Error: decayPion failed for mchi = $mchi"
#         continue
#     fi

#    # Check if the efficiency output file exists and is not empty
#           if [ -s "${efficiency_output_three_body[$k]}" ]; then
#               echo "Appending efficiency to ${total_efficiency_output_three_body[$k]}"
#               # Append the efficiency to the total output file
#               cat "${efficiency_output_three_body[$k]}" >> "${total_efficiency_output_three_body[$k]}"
#           else
#               echo "Warning: ${efficiency_output_three_body[$k]} is empty or does not exist. Skipping."
#           fi

#       done < mship_values.txt

#   done

# echo "All simulations completed. Total efficiencies are stored in ${total_efficiency_output_three_body[@]}"