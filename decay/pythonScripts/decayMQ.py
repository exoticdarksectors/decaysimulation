# import random
# import ROOT as r
# from math import *
#
# def Do2BodyDecay(p4_mother, m1, m2, cosTheta=random.uniform(-1,1), phi=random.uniform(-pi,pi)):
#     # get four-momenta p1,p2 of 2 daughter particles in decay mother -> 1 + 2
#     # p4_mother is a TLorentzVector of the four momentum of the mother particle
#     # m1, m2 are masses of 2 daughers
#     # cosTheta is cos(theta) of daughter 1 in the rest frame of the mother,
#     # where the z axis is along the original direction of the mother
#     # phi is phi in this system (probably just random between -pi and pi)
#     #
#     # returns a 2-tuple of TLorentzVectors p1,p2
#
#     p4_mother = p4_mother.Clone() # make sure we don't modify original
#
#     direction = p4_mother.BoostVector().Unit()
#     # if mother p4 is along z direction already, handle specially
#     if direction.Px() == 0.0  and direction.Py() == 0.0:
#         axis = r.TVector3(1,0,0)
#         angle = pi if direction.Pz() < 0 else 0.0
#     else:
#         axis = direction.Cross(r.TVector3(0,0,1))
#         angle = acos(direction.Dot(r.TVector3(0,0,1)))
#
#     p4_orig = p4_mother.Clone()
#     p4_mother.Rotate(angle, axis)
#
#     # print "({0:7.2f}, {1:7.2f}, {2:7.2f}), ({3:7.2f}, {4:7.2f}, {5:7.2f}) {6:7.2f}".format(p4_orig.Px(),p4_orig.Py(),p4_orig.Pz(),p4_mother.Px(),p4_mother.Py(),p4_mother.Pz(), p4_mother.P())
#
#     boost = p4_mother.BoostVector()
#     p4_mother.Boost(-boost)
#
#     theta = acos(cosTheta)
#     dir_1 = r.TVector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cosTheta)
#     dir_2 = -dir_1
#
#     M = p4_mother.M()
#     E1 = (M**2 + m1**2 - m2**2) / (2*M)
#     E2 = M - E1
#     p1 = sqrt(E1**2 - m1**2)
#     p2 = sqrt(E2**2 - m2**2)
#
#     p4_1 = r.TLorentzVector()
#     p4_1.SetPxPyPzE(p1*dir_1.x(), p1*dir_1.y(), p1*dir_1.z(), E1)
#     p4_2 = r.TLorentzVector()
#     p4_2.SetPxPyPzE(p2*dir_2.x(), p2*dir_2.y(), p2*dir_2.z(), E2)
#
#     p4_mother.Boost(boost)
#     p4_1.Boost(boost)
#     p4_2.Boost(boost)
#
#     p4_mother.Rotate(-angle, axis)
#     p4_1.Rotate(-angle, axis)
#     p4_2.Rotate(-angle, axis)
#
#     return p4_1, p4_2
#
#
# def DoDalitz(p4_mother, me, mX, useVME=True):
#     # Dalitz decay of the form P -> e+e-X (can be electron, but doesn't have to be)
#     # returns a 3-tuple of four-vectors pX, pe+, pe-
#
#     mP = p4_mother.M()
#
#     if 2*me + mX > mP:
#         raise ValueError("Illegal decay! 2*me + mX > mP for the decay P -> e+e-X")
#
#     # pdf of q^2 = m(e+e-)^2
#     print(useVME)
#     if useVME:
#         pdf = r.TF1("logq2_pdf","((1+exp(x)/([0]*[1]))^2-([0]+[1])^2*exp(x)/([0]*[1])^2)^1.5 * (1+0.5*[2]^2/exp(x)) * sqrt(1-[2]^2/exp(x)) * ([3]^4+([3]*[4])^2)/(([3]^2-exp(x))^2+([3]*[4])^2)", log((2*me)**2), log((mP-mX)**2))
#         pdf.SetParameter(0, mP-mX)
#         pdf.SetParameter(1, mP+mX)
#         pdf.SetParameter(2, 2*me)
#         pdf.SetParameter(3, 0.7755)  # part of the form factor F(q^2) = 1 + 0.03*q^2/mP^2
#         pdf.SetParameter(4, 0.1462)  # part of the form factor F(q^2) = 1 + 0.03*q^2/mP^2
#     else:
#         pdf = r.TF1("logq2_pdf", "((1+exp(x)/([0]*[1]))^2-([0]+[1])^2*exp(x)/([0]*[1])^2)^1.5 * (1+0.5*[2]^2/exp(x)) * sqrt(1-[2]^2/exp(x)) * (1+[3]*exp(x)/[0]^2)^2", log((2*me)**2), log((mP-mX)**2))
#         pdf.SetParameter(0, mP-mX)
#         pdf.SetParameter(1, mP+mX)
#         pdf.SetParameter(2, me)
#         pdf.SetParameter(3, 0.03)  # part of the form factor F(q^2) = 1 + 0.03*q^2/mP^2
#     pdf.SetNpx(1000)
#
#     q2 = exp(pdf.GetRandom())
#
#     # do the P -> X gstar decay
#     pX, pgstar = Do2BodyDecay(p4_mother, mX, sqrt(q2), random.uniform(-1,1), random.uniform(-pi,pi))
#
#     # pdf of cos(theta) in the gstar -> e+e- decay
#     pdf = r.TF1("cosTheta_pdf", "1 + x^2 + [0]^2/[1]*(1-x^2)", -1, 1)
#     pdf.SetParameter(0, 2*me)
#     pdf.SetParameter(1, q2)
#     pdf.SetNpx(1000)
#
#     cosTheta = pdf.GetRandom()
#     pe1, pe2 = Do2BodyDecay(pgstar, me, me, cosTheta, random.uniform(-pi,pi))
#
#     return pX, pe1, pe2
#
#
#
# if __name__=="__main__":
#
#     # for i in range(1):
#     #     p = r.TLorentzVector()
#     #     p.SetPtEtaPhiM(random.uniform(5,30), random.uniform(-3,3), random.uniform(-pi,pi), 0.139)
#     #     Do2BodyDecay(p, 0, 0, 0, 0)
#
#     p4_pi = r.TLorentzVector()
#     p4_pi.SetPtEtaPhiM(0, 0, 0, 0.139)
#     pX, pe1, pe2 = DoDalitz(p4_pi, 0.000511, 0.05)
#     pX.Print()
#     pe1.Print()
#     pe2.Print()
#     (pX+pe1+pe2).Print()
import random
import ROOT as r
from math import *

def Do2BodyDecay(p4_mother, m1, m2, cosTheta=None, phi=None):
    if cosTheta is None:
        cosTheta = random.uniform(-1, 1)
    if phi is None:
        phi = random.uniform(-pi, pi)

    p4_mother = p4_mother.Clone()  # make sure we don't modify the original

    direction = p4_mother.BoostVector().Unit()
    if direction.Px() == 0.0 and direction.Py() == 0.0:
        axis = r.TVector3(1, 0, 0)
        angle = pi if direction.Pz() < 0 else 0.0
    else:
        axis = direction.Cross(r.TVector3(0, 0, 1))
        angle = acos(direction.Dot(r.TVector3(0, 0, 1)))

    p4_orig = p4_mother.Clone()
    p4_mother.Rotate(angle, axis)

    boost = p4_mother.BoostVector()
    p4_mother.Boost(-boost)

    theta = acos(max(-1, min(1, cosTheta)))
    dir_1 = r.TVector3(sin(theta) * cos(phi), sin(theta) * sin(phi), cosTheta)
    dir_2 = -dir_1

    M = p4_mother.M()
    E1 = (M**2 + m1**2 - m2**2) / (2 * M)
    E2 = M - E1
    p1 = sqrt(max(0, E1**2 - m1**2))
    p2 = sqrt(max(0, E2**2 - m2**2))

    p4_1 = r.TLorentzVector()
    p4_1.SetPxPyPzE(p1 * dir_1.x(), p1 * dir_1.y(), p1 * dir_1.z(), E1)
    p4_2 = r.TLorentzVector()
    p4_2.SetPxPyPzE(p2 * dir_2.x(), p2 * dir_2.y(), p2 * dir_2.z(), E2)

    p4_mother.Boost(boost)
    p4_1.Boost(boost)
    p4_2.Boost(boost)

    p4_mother.Rotate(-angle, axis)
    p4_1.Rotate(-angle, axis)
    p4_2.Rotate(-angle, axis)

    return p4_1, p4_2