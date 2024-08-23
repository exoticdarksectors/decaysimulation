#ifndef MQDECAY_H
#define MQDECAY_H

#include <Math/LorentzVector.h>
#include <Math/PtEtaPhiM4D.h>
#include <TLorentzVector.h>
#include <TVector3.h>

/////////////////////////////////////////////////////////////
// define dGamma/dq^2 for dalitz decays in various ways here
// NOTE: need to multiply by 2 if X=gamma
/////////////////////////////////////////////////////////////
// dGamma / dlog(q^2), using VDM model
//    [0] is maximum q (mP-mX)
//    [1] is mP+mX
//    [2] is minimum q (2*me)
//    [3] is "mrho" in form factor (0.7755)
//    [4] is "width rho" in form factor (0.1462)
//    Range is from log((2*me)^2) to log((mP-mX)^2)
#define DGDLOGQ2_VDM "1/(137.036*3*pi)*((1+exp(x)/([0]*[1]))^2-([0]+[1])^2*exp(x)/([0]*[1])^2)^1.5 * (1+0.5*[2]^2/exp(x)) * sqrt(1-[2]^2/exp(x)) * ([3]^4+([3]*[4])^2)/(([3]^2-exp(x))^2+([3]*[4])^2)"
// dGamma / dlog(q^2), not using VDM model
//    [0] is maximum q (mP-mX)
//    [1] is mP+mX
//    [2] is minimum q (2*me)
//    [3] is "a" in the form factor F(q2) = 1+a*q2 / mpi^2 (0.03)
//    [4] is "mpi" in the form factor (0.1350)
//    Range is from log((2*me)^2) to log((mP-mX)^2)
#define DGDLOGQ2_NONVDM "1/(137.036*3*pi)*((1+exp(x)/([0]*[1]))^2-([0]+[1])^2*exp(x)/([0]*[1])^2)^1.5 * (1+0.5*[2]^2/exp(x)) * sqrt(1-[2]^2/exp(x)) * (1+[3]*exp(x)/[4]^2)^2"
// dGamma / dq^2, using VDM
#define DGDQ2_VDM "1/(137.036*3*pi*x)*((1+x/([0]*[1]))^2-([0]+[1])^2*x/([0]*[1])^2)^1.5 * (1+0.5*[2]^2/x) * sqrt(1-[2]^2/x) * ([3]^4+([3]*[4])^2)/(([3]^2-x)^2+([3]*[4])^2)"
// dGamma / dq^2, not using VDM
#define DGDQ2_NONVDM "1/(137.036*3*pi*x)*((1+x/([0]*[1]))^2-([0]+[1])^2*x/([0]*[1])^2)^1.5 * (1+0.5*[2]^2/x) * sqrt(1-[2]^2/x) * (1+[3]*x/[4]^2)^2"


// used to solve cubic equation to generate cos(theta) values
double pdf_ct(double ct, double a);
double cdf_ct(double ct, double a);
double newton(double a, double R, double guess=0.0, double tol=1e-4);


std::pair<TLorentzVector,TLorentzVector>
Do2BodyDecay(TLorentzVector p4_mother, double m1, double m2, double cosTheta=-999, double phi=-999);

std::tuple<TLorentzVector,TLorentzVector,TLorentzVector>
DoDalitz(TLorentzVector p4_mother, double me, double mX, bool useVDM=true);

template <class T>
TLorentzVector LVtoTLV(ROOT::Math::LorentzVector<T> p){
    TLorentzVector t = TLorentzVector();
    t.SetXYZM(p.x(), p.y(), p.z(), p.M());
    return t;
}

template <class T>
ROOT::Math::LorentzVector<T> TLVtoLV(TLorentzVector p){    
    return ROOT::Math::LorentzVector<T>(p.Pt(), p.Eta(), p.Phi(), p.M());
}

// support also ROOT::Math::LorentzVector
template <class T>
std::pair<ROOT::Math::LorentzVector<T>, ROOT::Math::LorentzVector<T> >
Do2BodyDecay(ROOT::Math::LorentzVector<T> p4_mother, double m1, double m2, double cosTheta=-999, double phi=-999){
    TLorentzVector p1,p2;
    std::tie(p1,p2) = Do2BodyDecay(LVtoTLV(p4_mother), m1, m2, cosTheta, phi);
    ROOT::Math::LorentzVector<T> q1 = TLVtoLV<T>(p1);
    ROOT::Math::LorentzVector<T> q2 = TLVtoLV<T>(p2);
    // return std::pair<ROOT::Math::LorentzVector<T>, ROOT::Math::LorentzVector<T> >(q1, q2);
    return std::pair<ROOT::Math::LorentzVector<T>, ROOT::Math::LorentzVector<T> >(q1, q2);
}

template <class T>
std::tuple<ROOT::Math::LorentzVector<T>, ROOT::Math::LorentzVector<T>, ROOT::Math::LorentzVector<T> >
DoDalitz(ROOT::Math::LorentzVector<T> p4_mother, double me, double mX, bool useVDM=true){
    TLorentzVector p1,p2,p3;
    std::tie(p1,p2,p3) = DoDalitz(LVtoTLV(p4_mother), me, mX, useVDM);
    ROOT::Math::LorentzVector<T> q1 = TLVtoLV<T>(p1);
    ROOT::Math::LorentzVector<T> q2 = TLVtoLV<T>(p2);
    ROOT::Math::LorentzVector<T> q3 = TLVtoLV<T>(p3);
    return std::make_tuple(q1, q2, q3);
}

#endif
