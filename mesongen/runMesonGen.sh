#!/bin/bash

# Define the input and output file paths

declare -a meson_id
declare -a meson_events

meson_id=(
    "111"  # pi0
    "221"  # eta
    "443"  # J/psi
    # "553"  # upsilon
    "113"  # rho
    "223"  # omega
    "333"  # phi
)

meson_events=(
    "50000"
    "200000"
    "100000000"
    # "100000"
    "1000000"
    "1000000"
    "10000000"
)

# Iterate over the array indices
for i in "${!meson_id[@]}"; do
    echo "Executing: cmake-build-debug/mesonGen 1234 10 ${meson_events[$i]} ${meson_id[$i]}"
    ./cmake-build-debug/mesonGen 1234 10 "${meson_events[$i]}" "${meson_id[$i]}"

    if [ $? -eq 0 ]; then
        echo "mesonGen execution successful."
    else
        echo "Error in mesonGen execution." >&2
        exit 1
    fi

    # Remove old merged file if it exists
    if [ -f "output-data/mesons_seed1234_${meson_id[$i]}.root" ]; then
        echo "File output-data/mesons_seed1234_${meson_id[$i]}.root exists. Removing..."
        rm "output-data/mesons_seed1234_${meson_id[$i]}.root"
    fi

    echo "Merging ROOT files with hadd..."
    hadd output-data/mesons_seed1234_${meson_id[$i]}.root output-data/mesons_seed1234_t{0..9}.root

    if [ $? -eq 0 ]; then
        echo "hadd successful."
    else
        echo "Error during hadd." >&2
        exit 1
    fi

    # --- Check if the ROOT file id matches the intended one ---
    echo "Checking meson id inside ROOT file..."
    root -l -b -q <<EOF | grep -q "ID check passed"
{
    TFile f("output-data/mesons_seed1234_${meson_id[$i]}.root");
    TTree* t = (TTree*)f.Get("mesons");
    int id;
    t->SetBranchAddress("id", &id);
    t->GetEntry(0);
    if (id == ${meson_id[$i]}) {
        std::cout << "ID check passed" << std::endl;
    } else {
        std::cerr << "ERROR: ROOT file meson ID (" << id
                  << ") does not match expected (${meson_id[$i]})." << std::endl;
        gSystem->Exit(1);
    }
}
EOF

    if [ $? -eq 0 ]; then
        echo "Meson ID matches expected ${meson_id[$i]}."
    else
        echo "Meson ID check failed for ${meson_id[$i]}." >&2
        exit 1
    fi
done
