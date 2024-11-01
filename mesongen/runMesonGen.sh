#!/bin/bash

# Define the input and output file paths

declare -a meson_id
declare -a meson_events

meson_id=(
    "111"  # pi0
    "221"  # eta
    "443"  #J/psi
    # "553"  # upsilon
    "113"  # rho
    "223"  # omega
    "333"  # phi
)

meson_events=(
    "50000"
    "200000"
    "10000000"
    # "100000"
    "1000000"
    "1000000"
    "10000000"
)

# Iterate over the array indices
for i in "${!meson_id[@]}"; do
    # Run the C++ script with the current mchi value
    echo "Executing: cmake-build-debug/mesonGen 1234 10 ${meson_events[$i]} ${meson_id[$i]}"
    ./cmake-build-debug/mesonGen 1234 10 "${meson_events[$i]}" "${meson_id[$i]}"

    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "mesonGen execution successful."
    else
        echo "Error in mesonGen execution." >&2
        exit 1
    fi

    # Check if the ROOT file exists and remove it if it does
    if [ -f "output-data/mesons_seed1234_${meson_id[$i]}.root" ]; then
        echo "File output-data/mesons_seed1234_${meson_id[$i]}.root exists. Removing..."
        rm "output-data/mesons_seed1234_${meson_id[$i]}.root"
    fi

    # Run hadd to merge the ROOT files
    echo "Merging ROOT files with hadd..."
    hadd output-data/mesons_seed1234_${meson_id[$i]}.root output-data/mesons_seed1234_t{0..9}.root

    # Check if hadd was successful
    if [ $? -eq 0 ]; then
        echo "hadd successful."
    else
        echo "Error during hadd." >&2
        exit 1
    fi
done
