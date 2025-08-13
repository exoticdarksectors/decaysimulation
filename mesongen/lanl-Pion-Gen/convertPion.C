// save this as convertPion.C
#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <iostream>

void convertPion() {
    // 1) Open output ROOT file
    TFile out("pion_production.root", "RECREATE");
    
    // 2) Create TTree and branches
    TTree tree("pion","pion kinematics");
    Float_t P, cosTheta, phi;
    tree.Branch("P",        &P,         "P/F");
    tree.Branch("cosTheta", &cosTheta,  "cosTheta/F");
    tree.Branch("phi",      &phi,       "phi/F");
    
    // 3) Open your ASCII file
    std::ifstream in("Pion_Production.txt");
    if (!in.is_open()) {
        std::cerr << "Error: cannot open Pion_Production.txt\n";
        return;
    }
    
    // 4) Loop over lines and fill tree
    while (in >> P >> cosTheta >> phi) {
        tree.Fill();
    }
    in.close();
    
    // 5) Write & close
    tree.Write();
    out.Close();
    
    std::cout << "Wrote " << tree.GetEntries() 
              << " entries into pion_production.root\n";
}