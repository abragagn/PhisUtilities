#ifndef PhisUtil_H
#define PhisUtil_H

#include <vector>
#include <string>
#include "PDAnalysis/Ntu/interface/PDEnumString.h"
#include "PDAnalysis/Ntu/interface/PDGenHandler.h"
#include "PDAnalyzerUtil.h"
#include "TF1.h"
#include "Math/Vector4D.h"
#include "Math/Vector3D.h"
#include "Math/Point3D.h"
#include "Math/MinimizerOptions.h"
#include "TMatrixD.h"
#include "TVectorD.h"

using namespace std;


class PhisUtil:  public virtual PDAnalyzerUtil
,                public virtual PDGenHandler{

public:

    PhisUtil();
    virtual ~PhisUtil();

    const double MUMASS  = 0.1056583745;
    const double EMASS   = 0.00051099895;
    const double PMASS   = 0.93827208816;

    const double PIMASS  = 0.13957039;
    const double PI0MASS = 0.1349768;

    const double KMASS   = 0.493677;
    const double K0MASS  = 0.497611;
    const double KXMASS  = 0.89166;

    const double BSMASS  = 5.36688;
    const double B0MASS  = 5.27965;
    const double BUMASS  = 5.27934;
    const double BCMASS  = 6.2749;

    const double JPSIMASS = 3.0969;
    const double PHIMASS  = 1.019461;

    const double LAMBDAMASS = 5.61960;

    const unsigned int listLundBmesons[24]  = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 
                                            20523, 515, 525, 531, 10531, 533, 10533, 20533, 535, 541, 
                                            10541, 543, 10543, 20543, 545};

    const unsigned int listLundBbaryons[35] = {5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 
                                            5312, 5322, 5314, 5324, 5332, 5334, 5142, 5242, 5412, 
                                            5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 
                                            5522, 5514, 5524, 5532, 5534, 5542, 5544, 5554};


    const unsigned int listLundBottonium[29] = {551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 
                                            20553, 30553, 100553, 110553, 120553, 130553, 200553, 210553, 
                                            220553, 300553, 9000553, 9010553, 555, 10555, 20555, 100555, 
                                            110555, 120555, 200555, 557, 100557};

    const unsigned int LongLivedList[7]     = {11,13,211,321,2212}; //e, mu, pi, K, p   

    const unsigned int listLundCmesons[18]  = {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 
                                            415, 425, 431, 10431, 433, 10433, 20433, 435};

    const unsigned int listLundCbaryons[22] = {4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 
                                            4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444};

    const unsigned int listLundCharmonium[13] = {441, 10441, 100441, 443, 10443, 20443, 100443, 30443, 9000443,
                                            9010443, 9020443, 445, 100445};

    double BsMassRange[2] = {5.0, 6.0};
    double BuMassRange[2] = {5.0, 6.0};
    double BdMassRange[2] = {5.0, 6.0};

    double BsMassRangeTight[2] = {5.24, 5.49};
    double BuMassRangeTight[2] = {5.00, 5.65};
    double BdMassRangeTight[2] = {5.00, 5.50};

    double bMassMin     = 5.24;
    double bMassMax     = 5.49;
    double jpsiMassWin  = 0.150;
    double phiMassWin   = 0.010;
    double kstarMassWin = 0.09;

    double bsPtMin      = 9.5;
    double muPtMin      = 3.5;  // change according to the hlt
    double kaonPtMin    = 1.15; // change according to the hlt
    double muEtaMax     = 2.4;
    double kaonEtaMax   = 2.5;
    int    kaonHitsMin  = 4;
    double vtxProbMin   = 0.02;
    int    softMuonID   = 1;

    double ctMin        = -999.;    // change according to the hlt
    double ctSigmaMin   = -999.;    // change according to the hlt

    const double muPtMin_jpsimu        = 3.5;
    const double muPtMin_jpsitrktrk    = 4.0;
    const double kaonPtMin_jpsimu      = 1.15;
    const double kaonPtMin_jpsitrktrk  = 0.9;    
    const double ctMin_jpsimu          = 0.006;
    const double ctMin_jpsitrktrk      = 0.01;
    const double ctSigmaMin_jpsimu     = -999.;
    const double ctSigmaMin_jpsitrktrk = 3.;


    // HLT
    double jpsiPt_HLTtk     = 6.9;
    double muPt_HLTtk       = 4.0;
    double kaonPt_HLTtk     = 0.8;
    double ctSigmaMin_HLTtk = 3.0;

    double muPt_HLTmu       = 3.5;


    void SetBMassMin(double newValue){ bMassMin = newValue; }
    void SetBMassMax(double newValue){ bMassMax = newValue; }
    void SetJpsiMassWin(double newValue){ jpsiMassWin = newValue; }
    void SetPhiMassWin (double newValue){ phiMassWin  = newValue; }
    void SetBsPtMin(double newValue){ bsPtMin = newValue; }
    void SetMuPtMin(double newValue){ muPtMin = newValue; }
    void SetKaonPtMin(double newValue){ kaonPtMin = newValue; }
    void SetMuEtaMax(double newValue){ muEtaMax = newValue; }
    void SetKaonEtaMax(double newValue){ kaonEtaMax = newValue; }
    void SetKaonHitsMin(double newValue){ kaonHitsMin = newValue; }
    void SetCtMin(double newValue){ ctMin = newValue; }
    void SetCtSigmaMin(double newValue){ ctSigmaMin = newValue; }
    void SetVtxProbMin(double newValue){ vtxProbMin = newValue; }
    void SetSoftMuonID(double newValue){ softMuonID = newValue; }

    void SetJpsiMuCut(bool ctCut);
    void SetJpsiTrkTrkCut(bool ctCut);

    void SetBsMassRange(double lower, double upper) { BsMassRange[0] = lower; BsMassRange[1] = upper; }
    void SetBuMassRange(double lower, double upper) { BuMassRange[0] = lower; BuMassRange[1] = upper; }
    void SetBdMassRange(double lower, double upper) { BdMassRange[0] = lower; BdMassRange[1] = upper; }

    bool IsLongLived( unsigned int i );
    bool IsB( unsigned int i ) ;
    bool IsC( unsigned int i ) ;
    bool IsCharmonium( unsigned int i ) ;
    bool IsBottomium( unsigned int i ) ;

    int    GetClosestGen( double eta, double phi, double pt );
    int    GetClosestGenLongLivedB( double eta, double phi, double pt, std::vector <int> *GenList );
    int    GetOverlappedTrack( int trk, std::vector <int> *List );
    bool   AreOverlapped( double pt1, double eta1, double phi1, double pt2, double eta2, double phi2 );
    int    GetAncestor( unsigned int iGen, std::vector <int> *GenList );
    int    MuonFromTrack(int trk);
    double GetGenCT( unsigned int genIndex );

    int    GetBestBstrange();
    int    GetBestBdown();
    int    GetBestBup();
    int    GetBestBstrangeTight();
    int    GetBestBdownTight();
    int    GetBestBupTight();
    bool   IsTightJPsi(int iJPsi);
    bool   IsTightPhi(int iPhi);
    int    GetBestJpsi();
    int    GetTightCandidate(TString process);
    int    GetCandidate(TString process);
    
    double GetInvMass(int i1, int i2, double mass1, double mass2);
    int    GetMixStatus( unsigned int genIndex );
    double GetMuoPFiso (int iMuon);
    double GetJetCharge(int iJet, double kappa);
    double GetListCharge(std::vector <int> *list, double kappa);
    int    IPsign(int iMuon, int iPV);
    double GetJetProbb(int iJet);
    double CountEventsWithFit(TH1 *hist, TString process);
    template <typename Vector>
    int    GetPVPointing(int isvt, Vector t);
    double dZ(int itk, int iPV);
    double dXYjet(int itk, int iPV, int iJet);
    void   PrintMotherChain(int iGen, std::ostream& out = std::cout);
    void   PrintDaughterTree(int iGen, const std::string & pre); // old interface
    void   PrintDaughterTreePt(int iGen, const std::string & pre); // old interface
    void   PrintDaughterTree(int iGen, std::ostream& out = std::cout, const std::string& prefix = "", bool isLastSibling = true);
    bool   HasDaughter(int iGen);

    template <typename Point, typename Vector>
    ROOT::Math::XYZPoint PCAwrtBeamSpot(Point point, Vector vector);

    template <typename Vector>
    double GetCt2D(Vector ptB, double bMass, int iSV);
    template <typename Vector>
    double GetCt2D(Vector ptB, double bMass, int iSV, int iPV, bool useRefittedPV);
    
    template <typename Vector>
    double GetCt3D(Vector ptB, double bMass, int iSV, int iPV, bool useRefittedPV);
    
    template <typename PtVector, typename LVector>
    double GetCtFromVector(PtVector ptB, double bMass, LVector L);
    
    template <typename Vector>
    double GetCt2DErr(Vector ptB, double bMass, int iSV);
    template <typename Vector>
    double GetCt2DErr(Vector ptB, double bMass, int iSV, int iPV, bool useRefittedPV);
    
    template <typename PtVector, typename LVector>
    double GetCtErrFromVector(PtVector ptB, double bMass, LVector L, TMatrixD PVSVcov, int iSV);
    // double GetCt3DErr(Vector ptB, int iSV, int iPV);

    bool isTrkHighPurity(int itk){ return (( trkQuality->at( itk ) >> 2 ) & 1); }
    int  HLTMatch(int, bool);
    
    template <typename Vector1, typename Vector2>
    double DeltaR(Vector1 v1, Vector2 v2);
    
    
    ROOT::Math::PxPyPzEVector svtP4( int iSV, bool fromRefit = true);


protected:


};

// =================================================================================================
// Finds the PCA between the beamspot line and the line passing from the SV parallel to the momentum
// Evaluated on the beamspot line
template <typename Point, typename Vector>
ROOT::Math::XYZPoint PhisUtil::PCAwrtBeamSpot(Point point, Vector vector)
{
    using namespace ROOT::Math;
    XYZPoint  sv(point);
    XYZVector p(vector);
    
    XYZPoint  bs(bsX, bsY, bsZ);
    XYZVector bsDir(bsdXdZ, bsdYdZ, 1);
    
    auto d = bs - sv;
    
    double numerator = d.Dot(p.Cross(p.Cross(bsDir)));
    
    double denominator = (bsDir.Mag2()*p.Mag2() - pow(p.Dot(bsDir), 2));
    
    double t = numerator/denominator;

    return bs + t * bsDir;
}

// =================================================================================================
template <typename PtVector, typename LVector>
double PhisUtil::GetCtFromVector(PtVector pB, double bMass, LVector L)
{
    ROOT::Math::XYZVector pBCart(pB);
    ROOT::Math::XYZVector LCart(L);
    
    return bMass*LCart.Dot(pBCart)/pBCart.Mag2();
}

// =================================================================================================
template <typename Vector>
double PhisUtil::GetCt2D(Vector t, double bMass, int iSV) 
{    
    using namespace ROOT::Math;
    XYZPoint SVpos( svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV) );

    // beam spot position need to be corrected as a function of z
    double Zpos = PCAwrtBeamSpot(SVpos, t).Z();
    double bsX_fix = bsX + bsdXdZ*(Zpos - bsZ);
    double bsY_fix = bsY + bsdYdZ*(Zpos - bsZ);

    XYZPoint PVpos(bsX_fix, bsY_fix, 0.); //bsZ makes no sense
    auto Lxy = SVpos - PVpos;
    Lxy.SetZ(0.);

    XYZVector ptB(t);
    ptB.SetZ(0.);

    return GetCtFromVector(ptB, bMass, Lxy);
}

// =================================================================================================
template <typename Vector>
double PhisUtil::GetCt2D(Vector t, double bMass, int iSV, int iPV, bool useRefittedPV)
{
    using namespace ROOT::Math;
    XYZPoint SVpos(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));

    XYZPoint PVpos;
    if(useRefittedPV) PVpos = XYZPoint(svtX->at(iPV),svtY->at(iPV),svtZ->at(iPV));
    else              PVpos = XYZPoint(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
    
    auto Lxy = SVpos - PVpos;
    Lxy.SetZ(0.);

    XYZVector ptB(t);
    ptB.SetZ(0.);
    
    return GetCtFromVector(ptB, bMass, Lxy);
}

// =================================================================================================
template <typename Vector>
double PhisUtil::GetCt3D(Vector t, double bMass, int iSV, int iPV, bool useRefittedPV)
{
    using namespace ROOT::Math;
    XYZPoint SVpos(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));

    XYZPoint PVpos;
    if(useRefittedPV) PVpos = XYZPoint(svtX->at(iPV),svtY->at(iPV),svtZ->at(iPV));
    else              PVpos = XYZPoint(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
    
    auto Lxyz = SVpos - PVpos;
    
    return GetCtFromVector(t, bMass, Lxyz);
}

// =====================================================================================
template <typename PtVector, typename LVector>
double PhisUtil::GetCtErrFromVector(PtVector ptB, double bMass, LVector L, TMatrixD PVSVcov, int iSV)
{
    using namespace ROOT::Math;

    double ct = GetCtFromVector(ptB, bMass, L);

    double ptBArray[] = {ptB.X(), ptB.Y(), ptB.Z()};
    TVectorD ptBVector(3, ptBArray);

    if(L.Mag2() == 0){
        cout << "Secondary vertex is exactly the same as PV" << endl;
        return -999.;
    }

    //Lxy contribution
    double ctErr_L2 = pow(bMass/ptB.Mag2(),2)*PVSVcov.Similarity(ptBVector);

    //Pt contribution
    auto& tkB  = tracksFromSV(iSV);
    int ntk = tkB.size();

    unordered_map<int,int> trkTotpp;

    for(int itpp = 0; itpp < nTrkPer; itpp++)
        trkTotpp[tppTrk->at(itpp)] = itpp;
    
    vector<XYZVector> trkMomentum(ntk);
    vector<int> qArray(ntk);
    vector<int> trk_idx(ntk);

    for( int i = 0; i < ntk; ++i ){
        int j = tkB.at(i);
        XYZVector tempTV3(trkPx->at(j), trkPy->at(j), trkPz->at(j));

        trkMomentum[i]  = tempTV3;
        qArray[i]       = trkCharge->at(j);
        trk_idx[i]      = j;
    }

    double ctErr_pt2 = 0.;

    // Construct trkTotpp map
    for(int i = 0; i < ntk; ++i){
        int itpp = trkTotpp.at(trk_idx[i]);

        double covHelixArray[] = {
            tppSQopQop->at(itpp), tppSQopLam->at(itpp), tppSQopPhi->at(itpp),
            tppSQopLam->at(itpp), tppSLamLam->at(itpp), tppSLamPhi->at(itpp),
            tppSQopPhi->at(itpp), tppSLamPhi->at(itpp), tppSPhiPhi->at(itpp)
        };
        TMatrixD covHelix(3,3);
        covHelix.SetMatrixArray(covHelixArray);

        // Helix parameters
        double Qop = qArray[i]/trkMomentum[i].R();
        double Lam = TMath::Pi()/2 - trkMomentum[i].Theta();

        double DQop = 1./(ptB.Mag2()*Qop)*(2*ct*ptB.Dot(trkMomentum[i]) - bMass*L.Dot(trkMomentum[i]));
        double DLam = tan(Lam)/ptB.Mag2()*(2*ct*ptB.Dot(trkMomentum[i]) - bMass*L.Dot(trkMomentum[i]));
        double DPhi = 1./ptB.Mag2()*(bMass*(L.Y()*trkMomentum[i].X() - L.X()*trkMomentum[i].Y())
                                  - 2*ct*(ptB.Y()*trkMomentum[i].X() - ptB.X()*trkMomentum[i].Y()));

        double DHelixArray[] = {DQop, DLam, DPhi};
        TVectorD DHelixVectorD(3, DHelixArray);

        ctErr_pt2 += covHelix.Similarity(DHelixVectorD);
    }
    
    return sqrt(ctErr_L2 + ctErr_pt2);
}

// =================================================================================================
template <typename Vector>
double PhisUtil::GetCt2DErr(Vector t, double bMass, int iSV, int iPV, bool useRefittedPV)
{
    using namespace ROOT::Math;
    XYZPoint SVpos(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));

    XYZPoint PVpos;
    if(useRefittedPV) PVpos = XYZPoint(svtX->at(iPV),svtY->at(iPV),svtZ->at(iPV));
    else              PVpos = XYZPoint(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
    
    auto Lxy = SVpos - PVpos;
    Lxy.SetZ(0.);

    XYZVector ptB(t);
    ptB.SetZ(0.);

    // Conversions for matrix multiplications
    double covSVArray[] = {
        svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
        svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV),
        svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)
    };

    array<double, 9> covPVArray;
    if(useRefittedPV){
        covPVArray = {{
            svtSxx->at(iPV),svtSxy->at(iPV),svtSxz->at(iPV),
            svtSxy->at(iPV),svtSyy->at(iPV),svtSyz->at(iPV),
            svtSxz->at(iPV),svtSyz->at(iPV),svtSzz->at(iPV)
        }};
    }else{
        covPVArray = {{
            pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
            pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV),
            pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)
        }};
    }

    TMatrixD covSV(3,3);
    TMatrixD covPV(3,3);
    covSV.SetMatrixArray(covSVArray);
    covPV.SetMatrixArray(covPVArray.data());

    TMatrixD covTot = covSV + covPV;
    
    return GetCtErrFromVector(ptB, bMass, Lxy, covTot, iSV);
}


// =================================================================================================
template <typename Vector>
double PhisUtil::GetCt2DErr(Vector t, double bMass, int iSV)
{
    using namespace ROOT::Math;
    XYZPoint SVpos(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));
    // beam spot position need to be corrected as a function of z
    double Zpos = PCAwrtBeamSpot(SVpos, t).Z();
    double bsX_fix = bsX + bsdXdZ*(Zpos - bsZ);
    double bsY_fix = bsY + bsdYdZ*(Zpos - bsZ);

    XYZPoint PVpos(bsX_fix, bsY_fix, 0.);
    auto Lxy = SVpos - PVpos;
    Lxy.SetZ(0.);

    XYZVector ptB(t);
    ptB.SetZ(0.);

    // Conversions for matrix multiplications
    double covSVArray[] = {
        svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
        svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV),
        svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)
    };

    double covPVArray[] = {
        pow(bwX,2),0.,0.,
        0.,pow(bwY,2),0.,
        0.,0.,0.
    };

    TMatrixD covSV(3,3);
    TMatrixD covPV(3,3);
    covSV.SetMatrixArray(covSVArray);
    covPV.SetMatrixArray(covPVArray);

    TMatrixD covTot = covSV + covPV;
    
    return GetCtErrFromVector(ptB, bMass, Lxy, covTot, iSV);
}

// =====================================================================================
template <typename Vector>
int PhisUtil::GetPVPointing(int iSV, Vector t)
{
    using namespace ROOT::Math;
    int bestPV = -1;
    double bestAngle = -1; // actually the cosine

    XYZPoint vSVT( svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV) );
    XYZVector vB(t);

    for( int i=0; i<nPVertices; ++i ){
       if(abs(svtZ->at(iSV) - pvtZ->at(i)) > 5. ) continue;

       XYZPoint vPV(pvtX->at(i), pvtY->at(i), pvtZ->at(i) );
       auto vPointing = vSVT - vPV;
       
       // genVector has no Angle method, apparently
       double angle = vPointing.Unit().Dot(vB.Unit());
       

       if(angle > bestAngle ){
         bestAngle = angle;
         bestPV = i;
       }
    }
    return bestPV;
}

// ========================================================================================
template <typename Vector1, typename Vector2>
double PhisUtil::DeltaR(Vector1 v1, Vector2 v2)
{
    double deta = v1.Eta() - v2.Eta();
    double dphi = TMath::Abs(TMath::Pi() - TMath::Abs(TMath::Pi() - TMath::Abs(v1.Phi() - v2.Phi())));
    
    return std::hypot(deta, dphi);
}

#endif
