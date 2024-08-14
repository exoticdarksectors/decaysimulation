// decay pion to mcp with massless dark photon mediator
// read pion momentum from root files generated by mesonGen
// save momentum distribution of mcp into mcp-production.root
// Author : Insung Hwang (his5624@korea.ac.kr)
// Co-author: Leo Bailloeul (lbailloeul@ucdavis.edu)
#include <iostream>

#include "TRandom3.h"
#include "TLorentzVector.h"
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include "TVector3.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace ROOT::Math;

// constants
const double mpi = 0.1349768; // GeV, charged pion mass
const double mp = 938.272;
const double alpha = 0.0072973526; // 1/137
const double PI = 3.14159265358979;
const double conv2rad = PI / 180.0;

double envelope = 150; // enveloping function of ddBrPi2gxx for rejection sampling

// model parameter, like coupling constants
const double epsilon = 1.0;
const double mchi = 0.0100; // GeV, mcp mass
const double BrPi2gg = 1e9; // Br(pion -> gamam gamma)

// detector parameters
double detectorRadius = 2; //meters
double distanceToBox = 40.0; //meters

// store particle information
struct Particle {
    double mass;
    double flightTime;
    TLorentzVector momentum;
    TLorentzVector startPosition;
    TLorentzVector endPosition;
};

//  sqrt(Kallen function(x,y,z))/2x = lambda, it defines momentum in two-body decay
double lambda(double x, double y, double z) {
    return sqrt(pow(x, 4) + pow(y, 4) + pow(z, 4) - 2 * pow(x * y, 2) - 2 * pow(y * z, 2) - 2 * pow(z * x, 2)) / (
               2 * x);
}

// consider only massless dark vector mediator, so no on-shell contribution
double ddBrPi2gxx(double s, double theta) {
    return sin(theta) * pow(epsilon, 2) * alpha / (4 * PI * s) * pow(1 - s / pow(mpi, 2), 3) *
           sqrt(1 - 4 * pow(mchi, 2) / s) * (2 - (1 - 4 * pow(mchi, 2) / s) * pow(sin(theta), 2)) * BrPi2gg;
}


int main(int argc, char* argv[]) {

    // Check if the correct number of arguments are provided
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <file1>" << endl;
        return 1;
    }


    // Input: open the ROOT file which was generated by mesonGen. It is a tree containing px, py, pz, and the mass m
    TFile *myFile = TFile::Open(argv[1]);
    if (!myFile || myFile->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return -1;
    }

    // Define input tree and branches reader
    TTree *inputTree = (TTree *) myFile->Get("mesons");
    double mypx, mypy, mypz, mye;
    inputTree->SetBranchAddress("px", &mypx);
    inputTree->SetBranchAddress("py", &mypy);
    inputTree->SetBranchAddress("pz", &mypz);
    inputTree->SetBranchAddress("e", &mye);

    // Define variables to receive tree parsing
    double momP, momTheta, momPhi;

    // Output: define the output tree mcp and its branches
    TFile *output = new TFile(argv[2], "RECREATE");
    TTree *tree = new TTree("mcp", "mcp");
    TTree *filteredTree = new TTree("mcp-filtered", "mcp-filtered");
    double PX, PY, PZ, PP, M, PHI, THETA;
    tree->Branch("Px", &PX);
    tree->Branch("Py", &PY);
    tree->Branch("Pz", &PZ);
    tree->Branch("PP", &PP);
    tree->Branch("Phi", &PHI);
    tree->Branch("Theta", &THETA);
    tree->Branch("M", &M);

    double fPX, fPY, fPZ, fPP, fM, fPHI, fTHETA;
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
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> lorentzVector(mypx, mypy, mypz, mye);
        // Get the azimuthal angle (phi)
        momPhi = lorentzVector.phi();
        // Get the polar angle (theta)
        momTheta = lorentzVector.theta();
        // Get Magnitude
        momP = lorentzVector.mag();

        double s, theta;
        while (true) {
            s = rand.Uniform(4.0 * pow(mchi, 2), pow(mpi, 2));
            theta = rand.Uniform(0.0, PI);
            double y = ddBrPi2gxx(s, theta);
            // cout << y << endl;
            if (y >= rand.Uniform(0, 1) * envelope) {
                // renew envelope when y > envelope
                if (y > envelope)
                    envelope = y;
                break;
            }
        }

        // two-body decay off-shell photon V -> xxbar in the rest frame of V
        double P = lambda(sqrt(s), mchi, mchi);
        double phi = rand.Uniform(0, 2.0 * PI);
        Particle mcp1;
        mcp1.momentum.SetPxPyPzE(P * sin(theta) * cos(phi), P * sin(theta) * sin(phi), P * cos(theta),
                                 sqrt(P * P + mchi * mchi));
        // cout << mcp1.momentum.Px() << " " << mcp1.momentum.Py() << " " << mcp1.momentum.Pz() << endl;
        Particle mcp2;
        mcp2.momentum.SetPxPyPzE(-P * sin(theta) * cos(phi), -P * sin(theta) * sin(phi), -P * cos(theta),
                                 sqrt(P * P + mchi * mchi));

        // boost V from its rest frame along z direction
        double vP = lambda(mpi, sqrt(s), 0); // pion -> V + gamma, mass of gamma = 0
        double vBeta = vP / sqrt(vP * vP + s);
        // cout << "vBeta : " << vBeta << endl;
        mcp1.momentum.Boost(0, 0, vBeta);
        mcp2.momentum.Boost(0, 0, vBeta);
        // pion -> V + gamma is an isotropic two-body decay processs
        double vTheta = rand.Uniform(0.0, PI);
        double vPhi = rand.Uniform(0.0, 2 * PI);
        // rotate z axis of V back to pion rest frame
        mcp1.momentum.RotateZ(vTheta);
        mcp1.momentum.RotateY(vPhi);
        mcp2.momentum.RotateZ(vTheta);
        mcp2.momentum.RotateY(vPhi);
        // boost pion rest frame to lab frame
        double momE = sqrt(pow(momP, 2) + pow(mpi, 2));
        double momBx = momP * sin(momTheta) * cos(momPhi) / momE;
        double momBy = momP * sin(momTheta) * sin(momPhi) / momE;
        double momBz = momP * cos(momTheta) / momE;
        mcp1.momentum.Boost(momBx, momBy, momBz);
        mcp2.momentum.Boost(momBx, momBy, momBz);

        // store momentum vector values to fill tree
        PX = mcp1.momentum.Px();
        PY = mcp1.momentum.Py();
        PZ = mcp1.momentum.Pz();
        PP = mcp1.momentum.P();
        M = mcp1.momentum.M();
        PHI = mcp1.momentum.Phi();
        THETA = mcp1.momentum.Theta();

        // set minimum requried angles for interection with detector
        double thetaRequired = atan2(detectorRadius,distanceToBox);
        double thetaFinal = THETA;
        // filter mcps intersecting detector: stored into an additional filtered tree
        if (thetaFinal <= thetaRequired) {
            fPX = mcp1.momentum.Px();
            fPY = mcp1.momentum.Py();
            fPZ = mcp1.momentum.Pz();
            fPP = mcp1.momentum.P();
            fM = mcp1.momentum.M();
            fPHI = mcp1.momentum.Phi();
            fTHETA = mcp1.momentum.Theta();
            filteredTree->Fill();
        }
        tree->Fill();

        // same process done for mcp2
        PX = mcp2.momentum.Px();
        PY = mcp2.momentum.Py();
        PZ = mcp2.momentum.Pz();
        PP = mcp2.momentum.P();
        M = mcp2.momentum.M();
        PHI = mcp2.momentum.Phi();
        THETA = mcp2.momentum.Theta();
        thetaFinal = THETA;

        // building filtered tree from mcp2
        if (thetaFinal <= thetaRequired) {
            fPX = mcp2.momentum.Px();
            fPY = mcp2.momentum.Py();
            fPZ = mcp2.momentum.Pz();
            fPP = mcp2.momentum.P();
            fM = mcp2.momentum.M();
            fPHI = mcp2.momentum.Phi();
            fTHETA = mcp2.momentum.Theta();
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

    return 0;
}