// decayPionText.cc
// Author: Insung Hwang & Leo Bailloeul
// π → V* → χχ model, input from text file "Pion_Production.txt" with columns px, py, pz (MeV)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

#include "TRandom3.h"
#include "TLorentzVector.h"
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace ROOT::Math;

// two‐body breakup momentum via Källén function
double lambda(double M, double m1, double m2) {
    return sqrt(pow(M,4) + pow(m1,4) + pow(m2,4)
               -2*pow(M*m1,2) -2*pow(m1*m2,2) -2*pow(m2*M,2))
         /(2.0*M);
}

// off‐shell π→V*(χχ) differential branching
double ddBrPi2gxx(double s, double theta, double mchi) {
    const double mpi     = 0.1349768;      // GeV
    const double alpha   = 0.0072973526;   // 1/137
    const double PI      = 3.141592653589793;
    const double epsilon = 1.0;
    const double BrPi2gg = 1e9;            // π→γγ branching

    double xi = 1.0 - 4.0*mchi*mchi/s;
    return sin(theta)
         * epsilon*epsilon * alpha / (4.0*PI*s)
         * pow(1.0 - s/(mpi*mpi), 3)
         * sqrt(xi) * (2.0 - xi*pow(sin(theta),2))
         * BrPi2gg;
}

// enveloping function for rejection sampling (will auto‐update)
double envelope = 150.0;

// lightweight particle container
struct Particle {
    TLorentzVector momentum;
};

int main(int argc, char** argv) {
    if (argc!=5) {
        cerr<<"Usage: "<<argv[0]
            <<" <output_root> <mchi_GeV> <eff_output_txt> <pion_text_file>\n";
        return 1;
    }
    const string outRoot    = argv[1];
    const double mchi       = atof(argv[2]);
    const string outEffFile = argv[3];
    const string infile     = argv[4];

    // open pion text file: columns px, py, pz (MeV)
    ifstream in(infile);
    if (!in) {
        cerr<<"Failed to open pion data "<<infile<<"\n";
        return 2;
    }

    // prepare output ROOT file and trees
    TFile *fout        = TFile::Open(outRoot.c_str(),"RECREATE");
    TTree *treeAll     = new TTree("mcp_all",  "All MCPs");
    TTree *treeDet1    = new TTree("mcp_det1", "Hits at detector #1");
    TTree *treeDet2    = new TTree("mcp_det2", "Hits at detector #2");

    // common branches
    double PX,PY,PZ,PP,M,PHI,THETA;
    auto branchify = [&](TTree* t){
        t->Branch("Px",&PX);
        t->Branch("Py",&PY);
        t->Branch("Pz",&PZ);
        t->Branch("P" ,&PP);
        t->Branch("M" ,&M );
        t->Branch("Phi",&PHI);
        t->Branch("Theta",&THETA);
    };
    branchify(treeAll);
    branchify(treeDet1);
    branchify(treeDet2);

    // detector parameters (same radius, two distances)
    const double detectorRadius = 0.28209479176; // m
    const double distance1      = 6.0;           // m
    const double distance2      = 35.0;          // m
    const double thetaMax1      = atan2(detectorRadius, distance1);
    const double thetaMax2      = atan2(detectorRadius, distance2);

    TRandom3 rnd(0);

    // loop over pions: read px,py,pz in MeV
    double px_MeV, py_MeV, pz_MeV;
    double pVal, thVal, phVal;
    int count = 0;
    while (in >> px_MeV >> py_MeV >> pz_MeV) {
        // convert to spherical coords
        double p_lab = sqrt(px_MeV*px_MeV + py_MeV*py_MeV + pz_MeV*pz_MeV);
        if (p_lab <= 0) continue;
        pVal  = p_lab / 1000.0;              // GeV
        thVal = acos(pz_MeV / p_lab);        // [0, π]
        phVal = atan2(py_MeV, px_MeV);       // [–π, π]
        ++count;

        // sample off‐shell mass s and decay angle θ_rel
        double s, th_rel;
        while (true) {
            s      = rnd.Uniform(4*mchi*mchi, 0.1349768*0.1349768);
            th_rel = rnd.Uniform(0.0, M_PI);
            double y = ddBrPi2gxx(s, th_rel, mchi);
            if (y >= rnd.Uniform(0.0, envelope)) {
                if (y > envelope) envelope = y;
                break;
            }
        }

        // V*→χχ in V-rest frame
        double Pdec = lambda(sqrt(s), mchi, mchi);
        double phD  = rnd.Uniform(0, 2*M_PI);
        Particle a,b;
        a.momentum.SetPxPyPzE(
          Pdec*sin(th_rel)*cos(phD),
          Pdec*sin(th_rel)*sin(phD),
          Pdec*cos(th_rel),
          sqrt(Pdec*Pdec + mchi*mchi)
        );
        b.momentum = a.momentum;
        b.momentum.SetPx(-a.momentum.Px());
        b.momentum.SetPy(-a.momentum.Py());
        b.momentum.SetPz(-a.momentum.Pz());

        // isotropic π→V+γ boost back to π-rest
        double vP    = lambda(0.1349768, sqrt(s), 0.0);
        double vBeta = vP / sqrt(vP*vP + s);
        double vTh   = rnd.Uniform(0.0, M_PI);
        double vPh   = rnd.Uniform(0.0, 2*M_PI);
        a.momentum.Boost(
          vBeta*sin(vTh)*cos(vPh),
          vBeta*sin(vTh)*sin(vPh),
          vBeta*cos(vTh)
        );
        b.momentum = a.momentum; // same boost

        // boost π-rest → lab
        double Eπ = sqrt(pVal*pVal + 0.1349768*0.1349768);
        double bx = pVal*sin(thVal)*cos(phVal)/Eπ;
        double by = pVal*sin(thVal)*sin(phVal)/Eπ;
        double bz = pVal*cos(thVal)/Eπ;
        a.momentum.Boost(bx,by,bz);
        b.momentum.Boost(bx,by,bz);

        // fill both MCPs
        for (Particle* chi : { &a, &b }) {
            const TLorentzVector &mom = chi->momentum;
            PX    = mom.Px();
            PY    = mom.Py();
            PZ    = mom.Pz();
            PP    = mom.P();
            M     = mom.M();
            PHI   = mom.Phi();
            THETA = mom.Theta();

            treeAll->Fill();
            if (PZ > 0) {
                double x1 = distance1 * PX / PZ;
                double y1 = distance1 * PY / PZ;
                if (x1*x1 + y1*y1 <= detectorRadius*detectorRadius)
                    treeDet1->Fill();

                double x2 = distance2 * PX / PZ;
                double y2 = distance2 * PY / PZ;
                if (x2*x2 + y2*y2 <= detectorRadius*detectorRadius)
                    treeDet2->Fill();
            }
        }
    }
    in.close();
    cout<<"\nProcessed "<<count<<" pions.\n";

    // compute efficiencies
    double tot = treeAll->GetEntries();
    double n1  = treeDet1->GetEntries();
    double n2  = treeDet2->GetEntries();
    double eff1= n1/tot;
    double eff2= n2/tot;
    cout<<fixed<<setprecision(6)
        <<"Det1 eff = "<<eff1
        <<", Det2 eff = "<<eff2<<"\n";

    // write ROOT and text outputs
    fout->Write();
    fout->Close();

    ofstream foutTxt(outEffFile);
    foutTxt<<setprecision(6)
           <<"mchi "<<mchi
           <<"  eff1 "<<eff1 << n1
           <<"  eff2 "<<eff2 <<n2<<"\n";
    foutTxt.close();
    cout<<"Wrote efficiencies to "<<outEffFile<<"\n";

    return 0;
}
