// decay pion to mcp with massless dark photon mediator
// read pion momentum from root files generated by mesonGen
// save momentum distribution of mcp into mcp-production.root
// Author : Leo Bailloeul (lbailloeul@ucdavis.edu)
// Do2BodyDecay function taken from https://github.com/milliQan-sw/milliq_mcgen/blob/master/utils/decay.cc

#include <iostream>
#include "TRandom3.h"
#include "TLorentzVector.h"
#include <Math/GenVector/LorentzVector.h>
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include <TF1.h>
#include <TRandom.h>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;
using namespace ROOT::Math;

// constants
const double PI = 3.14159265358979;
const double conv2rad = PI / 180.0;
// model parameter, like coupling constants
// const double mchi = 0.0450;
// detector parameters
double detectorRadius = 1.0; //meters
double distanceToBox = 40.0; //meters

// use global TF1s. Otherwise, the overhead of re-computing integrals every time is verrrry slow
TF1 *PDF_LOGQ2_VDM = 0;
TF1 *PDF_LOGQ2_NONVDM = 0;

std::pair<TLorentzVector,TLorentzVector> Do2BodyDecay(TLorentzVector p4_mother, double m1, double m2, double cosTheta, double phi){
    // get four-momenta p1,p2 of 2 daughter particles in decay m -> d1 + d2
    // p4_mother is four momentum of mother particle
    // m1, m2 are masses of daughters d1 and d2
    // cosTheta is cos(theta) of d1 in the rest frame of the mother,
    // measured w.r.t. original direction of mother
    // phi is phi in this system.
    // If these are not provided, they are generated randomly in ranges (-1,1), (-pi,pi).
    // returns a pair of four-momenta p1,p2 of the daughters d1,d2

    if(m1+m2 > p4_mother.M()){
        std::cout << "ERROR: illegal 2-body decay! m1 + m2 > M (" << m1 << ", " << m2 << ", " << p4_mother.M() << ")\n";
        throw std::exception();
    }

    TVector3 direction = p4_mother.BoostVector().Unit();
    TVector3 axis;

    double angle;
    // special handling for the case where mother p4 is already along z-direction
    if(direction.Px()==0.0 && direction.Py()==0.0){
        axis = TVector3(1,0,0);
        angle = direction.Pz() < 0 ? M_PI : 0.0;
    }else{
        axis = direction.Cross(TVector3(0,0,1));
        angle = acos(direction.Dot(TVector3(0,0,1)));
    }

    // rotate mother so it points along +z axis
    p4_mother.Rotate(angle, axis);

    // boost mother so that it is at rest
    TVector3 boost = p4_mother.BoostVector();
    p4_mother.Boost(-boost);

    // assign cosTheta/phi randomly if they weren't provided
    if(cosTheta < -998)
        cosTheta = gRandom->Uniform(-1, 1);
    if(phi < -998)
        phi = gRandom->Uniform(-M_PI, M_PI);

    double theta = acos(cosTheta);
    TVector3 dir_1 = TVector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cosTheta);
    TVector3 dir_2 = -dir_1;

    double M = p4_mother.M();
    double E1 = (M*M + m1*m1 - m2*m2) / (2*M);
    double E2 = M - E1;
    double p1 = sqrt(E1*E1 - m1*m1);
    double p2 = sqrt(E2*E2 - m2*m2);

    TLorentzVector p4_1, p4_2;
    p4_1.SetPxPyPzE(p1*dir_1.x(), p1*dir_1.y(), p1*dir_1.z(), E1);
    p4_2.SetPxPyPzE(p2*dir_2.x(), p2*dir_2.y(), p2*dir_2.z(), E2);

    p4_1.Boost(boost);
    p4_2.Boost(boost);
    p4_1.Rotate(-angle, axis);
    p4_2.Rotate(-angle, axis);

    return std::pair<TLorentzVector,TLorentzVector> (p4_1, p4_2);
}

int main(int argc, char* argv[]) {

    // Check if the correct number of arguments are provided
    // input arguments for this scripts: ./decayVectorMeson.cc Meson_inputRootfile MCP_outputRootfile MCP_mass single_efficiency_value_outputfile
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " <file1>" << endl;
        return 1;
    }

    // Input: open the ROOT file which was generated by mesonGen. It is a tree containing px, py, pz, and the mass m
    TFile *myFile = TFile::Open(argv[1]);
    if (!myFile || myFile->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return -1;
    }
    double mchi = atof(argv[3]);

    // Define input tree and branches reader
    TTree *inputTree = (TTree *) myFile->Get("mesons");
    double mypx, mypy, mypz, mye;
    inputTree->SetBranchAddress("px", &mypx);
    inputTree->SetBranchAddress("py", &mypy);
    inputTree->SetBranchAddress("pz", &mypz);
    inputTree->SetBranchAddress("e", &mye);

    // Output: define the output tree mcp and its branches
    TFile *output = new TFile(argv[2], "RECREATE");
    TTree *tree = new TTree("mcp", "mcp");
    TTree *filteredTree = new TTree("mcp-filtered", "mcp-filtered");
    double PX, PY, PZ, PP, M, PHI, THETA, E;
    tree->Branch("Px", &PX);
    tree->Branch("Py", &PY);
    tree->Branch("Pz", &PZ);
    tree->Branch("PP", &PP);
    tree->Branch("Phi", &PHI);
    tree->Branch("Theta", &THETA);
    tree->Branch("M", &M);

    double fPX, fPY, fPZ, fPP, fM, fPHI, fTHETA, fE;
    filteredTree->Branch("Px", &fPX);
    filteredTree->Branch("Py", &fPY);
    filteredTree->Branch("Pz", &fPZ);
    filteredTree->Branch("PP", &fPP);
    filteredTree->Branch("Phi", &fPHI);
    filteredTree->Branch("Theta", &fTHETA);
    filteredTree->Branch("M", &fM);

    TRandom3 rand;

    // Loop through tree data
    int totalSize = inputTree->GetEntries();
    cout << "Number of input tree entries: " << totalSize << endl;

    for (int i = 0; i < totalSize; ++i) {
        // Display progress information
        if ((totalSize > 10000) && (i % (totalSize / 50) == 0)) {
            cout << "\rCompletion Percentage: " << fixed << setprecision(2) << 100.0 * (double) i / (double) totalSize
                    << "%" << flush;
        }

        // Get tree node entry
        inputTree->GetEntry(i);

        // Initialize Lorentz vector with vector components and particle energy
        TLorentzVector motherParticle(mypx, mypy, mypz, mye);

        double costheta = rand.Uniform(-1, 1);
        double phi = rand.Uniform(0, 2.0 * PI);

        std::pair<TLorentzVector, TLorentzVector> decayProducts = Do2BodyDecay(motherParticle, mchi, mchi, costheta, phi);

        // Extract the individual TLorentzVectors from the pair
        TLorentzVector mcp1 = decayProducts.first;
        TLorentzVector mcp2 = decayProducts.second;

        // store momentum vector values to fill tree
        PX = mcp1.Px();
        PY = mcp1.Py();
        PZ = mcp1.Pz();
        PP = mcp1.P();
        M = mcp1.M();
        PHI = mcp1.Phi();
        THETA = mcp1.Theta();

        // set minimum requried angles for interection with detector
        double thetaRequired = atan2(detectorRadius,distanceToBox);
        double thetaFinal = THETA;
        // filter mcps intersecting detector: stored into an additional filtered tree
        if (thetaFinal <= thetaRequired) {
            fPX = mcp1.Px();
            fPY = mcp1.Py();
            fPZ = mcp1.Pz();
            fPP = mcp1.P();
            fM = mcp1.M();
            fPHI = mcp1.Phi();
            fTHETA = mcp1.Theta();
            filteredTree->Fill();
        }
        tree->Fill();

        // same process done for mcp2
        PX = mcp2.Px();
        PY = mcp2.Py();
        PZ = mcp2.Pz();
        PP = mcp2.P();
        M = mcp2.M();
        PHI = mcp2.Phi();
        THETA = mcp2.Theta();
        thetaFinal = THETA;

        // building filtered tree from mcp2
        if (thetaFinal <= thetaRequired) {
            fPX = mcp2.Px();
            fPY = mcp2.Py();
            fPZ = mcp2.Pz();
            fPP = mcp2.P();
            fM = mcp2.M();
            fPHI = mcp2.Phi();
            fTHETA = mcp2.Theta();
            filteredTree->Fill();
        }
        tree->Fill();
    }

    // Move to the next line after completion
    cout << endl;

    // calculate efficiency: (filtered entry size)/(original entry size)
    int treeSize = tree->GetEntries();
    int filteredTreeSize = filteredTree->GetEntries();
    double efficiency = 100.0 * (double) filteredTreeSize/ (double) treeSize;
    cout << "Efficiency: " << efficiency << " %"<< endl;
    // Close the ROOT file
    myFile->Close();

    // Write in output ROOT file
    output->Write("", TObject::kOverwrite);
    output->Close();

    cout << "Completed Successfully!" << endl;
    cout << "Output stored in: " << argv[2] << endl;

    std::string output_filename = "../sensitivity-plot/" + std::string(argv[4]);

    // Create an ofstream object to write to the file
    std::ofstream output_file(output_filename);

    // Check if the file is open
    if (output_file.is_open()) {
        // Write the efficiency value to the file
        output_file << argv[3] << " " << std::to_string((efficiency/100.0)) << std::endl;

        // Close the file
        output_file.close();
        std::cout << "Efficiency written to " << output_filename << std::endl;
    } else {
        std::cerr << "Error: Unable to open file " << output_filename << std::endl;
        return 1;
    }

    return 0;
}
