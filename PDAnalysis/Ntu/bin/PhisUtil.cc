//////////////////////////////////////////////////////////////////////////////////
//                                                                              //
//      Collection of utility functions for the Bs->J/PsiPhi CPV analysis       //
//                                                                              //
//      Author: Alberto Bragagnolo (alberto.bragagnolo@cern.ch)                 //
//                                                                              //
//////////////////////////////////////////////////////////////////////////////////

#include "PhisUtil.h"

#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))

PhisUtil::PhisUtil() {}

PhisUtil::~PhisUtil() {}

using namespace std;

// =====================================================================================
bool PhisUtil::IsB( uint genindex ) 
{
    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundBmesons); ++i ) if( genCode == listLundBmesons[i] )   return true;
    for( uint i=0; i<ARRAY_SIZE(listLundBbaryons); ++i ) if( genCode == listLundBbaryons[i] ) return true;
    return false;
}

// =================================================================================================
bool PhisUtil::IsBottomium(uint genindex)
{
    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundBottonium); ++i ) if( genCode == listLundBottonium[i] ) return true;
    return false;
}

// =====================================================================================
bool PhisUtil::IsC( uint genindex ) 
{
    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundCmesons); ++i ) if( genCode == listLundCmesons[i] )   return true;
    for( uint i=0; i<ARRAY_SIZE(listLundCbaryons); ++i ) if( genCode == listLundCbaryons[i] ) return true;
    return false;
}

// =================================================================================================
bool PhisUtil::IsCharmonium(uint genindex)
{
    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(listLundCharmonium); ++i ) if( genCode == listLundCharmonium[i] ) return true;
    return false;
}

// =====================================================================================
bool PhisUtil::IsLongLived( uint genindex ) 
{
    uint genCode = abs( genId->at(genindex) );
    for( uint i=0; i<ARRAY_SIZE(LongLivedList); ++i ) if( genCode == LongLivedList[i] ) return true;
    return false;
}

// =====================================================================================
int PhisUtil::GetClosestGen( double eta, double phi, double pt ) 
{
    double drb = 0.12;
    double dpb = 0.3; 
    int best = -1;
    
    for( uint i=0; i<genId->size(); ++i ){
       if( !IsLongLived(i) ) continue;
       double dr = deltaR(eta, phi, genEta->at(i), genPhi->at(i));
       double dpt = abs(genPt->at(i) - pt)/genPt->at(i);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = (int) i;
       drb = dr;
    } 

    return best;
}

// =====================================================================================
int PhisUtil::GetOverlappedTrack( int trk, vector<int> *List )
{
    double drb = 0.1;
    double dpb = 0.1; 
    int best = -1;
    
    for(int it:*List){

       double dr = deltaR(trkEta->at(trk), trkPhi->at(trk), trkEta->at(it), trkPhi->at(it));
       double dpt = abs(trkPt->at(it) - trkPt->at(trk))/trkPt->at(it);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = it;
       drb = dr;
    }

    return best;
}

// =====================================================================================
bool PhisUtil::AreOverlapped( double pt1, double eta1, double phi1, double pt2, double eta2, double phi2 )
{
    double drb = 0.01;
    double dpb = 0.05; 

    double dr = deltaR(eta1, phi1, eta2, phi2);
    double dpt = abs(pt1 - pt2)/pt1;

    if( dr > drb ) return false;
    if( dpt > dpb) return false;

    return true;
}

// =====================================================================================
int PhisUtil::GetClosestGenLongLivedB( double eta, double phi, double pt, vector<int> *GenList ) 
{
    double drb = 0.4;
    double dpb = 0.4; 
    int best = -1;
    
    for(int it:*GenList){

       double dr = deltaR(eta, phi, genEta->at(it), genPhi->at(it));
       double dpt = abs(genPt->at(it) - pt)/genPt->at(it);

       if( dr > drb ) continue;
       if( dpt > dpb) continue;

       best = it;
       drb = dr;
    } 

    return best;
}

// =====================================================================================
int PhisUtil::GetAncestor( uint iGen, vector<int> *GenList ) 
{
    const vector<int>* aM = &allMothers(iGen);
    while( aM->size()>0 ){ 
       int a = aM->at(0);
       for( uint i=0; i<GenList->size(); ++i ) if( GenList->at(i) == a ) return GenList->at(i);
       aM = &allMothers( aM->at(0) );
    }
    return -1;
}

// =====================================================================================
int PhisUtil::MuonFromTrack(int trk)
{
    for( int iMuon = 0; iMuon<nMuons; ++iMuon){
       if(muonTrack( iMuon, PDEnumString::muInner ) == trk) return iMuon;
    }
    return -1;
}

// =====================================================================================
double PhisUtil::GetGenCT( uint genIndex ) 
{
    const vector<int>& aD = allDaughters(genIndex);
    if( aD.size() == 0 ) return -1;

    uint mthIndex = aD[0];

    if( genId->at( genIndex ) == - genId->at(genMother->at(genIndex)) ) mthIndex = genMother->at(genIndex ); 
 
    ROOT::Math::PtEtaPhiMVector pGen(genPt->at(genIndex), genEta->at(genIndex), genPhi->at(genIndex), genMass->at(genIndex));

    double dx = genVx->at(genIndex)-genVx->at(mthIndex);
    double dy = genVy->at(genIndex)-genVy->at(mthIndex);
    double dz = genVz->at(genIndex)-genVz->at(mthIndex);

    return sqrt( dx*dx+dy*dy+dz*dz )/pGen.Beta()/pGen.Gamma();
}

// ========================================================================================
int PhisUtil::GetBestBstrange()
{
    int index = -1;
    double bestVtxProb = 0.;
    for( int iB=0; iB<nSVertices; ++iB ){
       if((svtType->at(iB)!=PDEnumString::svtBsJPsiPhi) ) continue;
       if( svtMass->at(iB)<BsMassRange[0] || svtMass->at(iB)>BsMassRange[1] ) continue;

       double vtxprob = TMath::Prob(svtChi2->at(iB),svtNDOF->at(iB));
       if( vtxprob < bestVtxProb ) continue;

       index = iB;
       bestVtxProb = vtxprob;
    }
    return index;
}

// ========================================================================================
int PhisUtil::GetBestBup()
{
    int index = -1;
    double bestVtxProb = 0.;
    for( int iB=0; iB<nSVertices; ++iB ){
        if((svtType->at(iB)!=PDEnumString::svtBuJPsiK) ) continue;
        if( svtMass->at(iB)<BuMassRange[0] || svtMass->at(iB)>BuMassRange[1] ) continue;

        double vtxprob = TMath::Prob(svtChi2->at(iB),svtNDOF->at(iB));
        if( vtxprob < bestVtxProb ) continue;

        index = iB;
        bestVtxProb = vtxprob;
    }
    return index;
}

// ========================================================================================
int PhisUtil::GetBestBdown()
{
    int index = -1;
    double bestVtxProb = 0.;
    for( int iB=0; iB<nSVertices; ++iB ){
       if((svtType->at(iB)!=PDEnumString::svtBdJPsiKx) ) continue;
       if( svtMass->at(iB)<BdMassRange[0] || svtMass->at(iB)>BdMassRange[1] ) continue;

       double vtxprob = TMath::Prob(svtChi2->at(iB),svtNDOF->at(iB));
        if( vtxprob < bestVtxProb ) continue;

        index = iB;
        bestVtxProb = vtxprob;
    }
    return index;
}

// ========================================================================================
bool PhisUtil::IsTightJPsi(int iJPsi)
{
    if(abs(svtMass->at(iJPsi) - JPSIMASS) > jpsiMassWin ) return false;

    vector<int> tkJpsi = tracksFromSV(iJPsi);
    // ROOT::Math::PxPyPzEVector tJPsi(0,0,0,0);           // No more requiremebt on Jpsi momentum

    for( uint i=0; i<tkJpsi.size(); ++i ){
        int j = tkJpsi[i];
        if(trkPt->at(j) < muPtMin) return false;
        if(abs(trkEta->at(j)) > muEtaMax) return false;
        // ROOT::Math::PtEtaPhiMVector a( trkPt->at(j), trkEta->at(j), trkPhi->at(j), MUMASS );
        // tJPsi += a;
    }

    //if(tJPsi.Pt() < 7.0) return false;

    return true;
}

// ========================================================================================
bool PhisUtil::IsTightPhi(int iPhi)
{
    if(abs(svtMass->at(iPhi) - PHIMASS) > phiMassWin ) return false;
    
    vector<int> tkPhi = tracksFromSV(iPhi);
    for( uint i=0; i<tkPhi.size(); ++i ){
        int j = tkPhi[i];
        if(trkPt->at(j) < kaonPtMin) return false;
        if(abs(trkEta->at(j)) > kaonEtaMax) return false;

        int K_Hits = trkHitPattern->at(tkPhi[i]);
        if(((int(K_Hits)/100)%10000)%100<kaonHitsMin ) return false;
    }
    return true;
}

// ========================================================================================
int PhisUtil::GetBestBstrangeTight()
{
    int index = -1;
    double best = 0.;

    for( int iB = 0; iB < nSVertices; ++iB ){

        if((svtType->at(iB)!=PDEnumString::svtBsJPsiPhi) ) continue;

        int iJPsi = (daugSVsFromVtx(iB, PDEnumString::vtrCascade)).at(0);
        int iPhi  = (daugSVsFromVtx(iB, PDEnumString::vtrCascade)).at(1);

        vector<int> tkSsB  = tracksFromSV(iB);
        vector<int> tkJpsi = tracksFromSV(iJPsi);
        vector<int> tkPhi  = tracksFromSV(iPhi);

        //JPSI
        if(!IsTightJPsi(iJPsi)) continue;

        //PHI
        if(!IsTightPhi(iPhi)) continue;

        //BS
        ROOT::Math::PxPyPzEVector tB(0,0,0,0);

        for( uint iSsBTrk = 0; iSsBTrk < tkSsB.size(); ++iSsBTrk ){
            int j = tkSsB[iSsBTrk];
            double m = KMASS;
            if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MUMASS;
            ROOT::Math::PtEtaPhiMVector a( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
        }

        double bVprob = TMath::Prob( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BsMassRangeTight[0] || svtMass->at(iB)>BsMassRangeTight[1] ) continue;
        if( bVprob < vtxProbMin ) continue;
        if(tB.Pt() < bsPtMin) continue;

        // int iPV = GetPVPointing(iB, tB);
        // if(iPV<0) continue;
        if(GetCt2D(tB, BSMASS, iB) < ctMin) continue;
        if(GetCt2D(tB, BSMASS, iB) / GetCt2DErr(tB, BSMASS, iB) < ctSigmaMin) continue;

        if( bVprob < best ) continue;
        index = iB;
        best = bVprob;

    }
    return index;
}

// ========================================================================================
int PhisUtil::GetBestBupTight()
{
    int index = -1;
    double best = 0.;
    for( int iB=0; iB<nSVertices; ++iB ){

        if((svtType->at(iB)!=PDEnumString::svtBuJPsiK) ) continue;

        int iJPsi = (daugSVsFromVtx(iB, PDEnumString::vtrCascade)).at(0);
        if(!IsTightJPsi(iJPsi)) continue;

        vector<int> tkJpsi = tracksFromSV(iJPsi);
        vector<int> tkSsB = tracksFromSV(iB);

        ROOT::Math::PxPyPzEVector tB(0,0,0,0);
        double KaonPt = 0.;

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            double m = KMASS;
            if( j == tkJpsi[0] || j == tkJpsi[1] ){ m = MUMASS; }else{ KaonPt = trkPt->at(j); }
            ROOT::Math::PtEtaPhiMVector a( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
       }

       //Bu
        double bVprob = TMath::Prob( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BuMassRangeTight[0] || svtMass->at(iB)>BuMassRangeTight[1] ) continue;
        if( bVprob < vtxProbMin ) continue;
        if(tB.Pt() < bsPtMin) continue;
        if(KaonPt < 1.6) continue;

        // int iPV = GetPVPointing(iB, tB);
        // if(iPV<0) continue;
        if(GetCt2D(tB, BUMASS, iB) < ctMin) continue;
        if(GetCt2D(tB, BUMASS, iB) / GetCt2DErr(tB, BUMASS, iB) < ctSigmaMin) continue;

        if( bVprob < best ) continue;
        index = iB;
        best = bVprob;
    }
    return index;
}

// ========================================================================================
int PhisUtil::GetBestBdownTight()
{
    int index = -1;
    double best = 0.;
    for( int iB=0; iB<nSVertices; ++iB ){

        if((svtType->at(iB)!=PDEnumString::svtBdJPsiKx) ) continue;

        int iJPsi = (daugSVsFromVtx(iB, PDEnumString::vtrCascade)).at(0);
        if(!IsTightJPsi(iJPsi)) continue;

        vector<int> tkJpsi = tracksFromSV(iJPsi);
        vector<int> tkSsB = tracksFromSV(iB);

        ROOT::Math::PxPyPzEVector tB(0,0,0,0);

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            double m = KXMASS;
            if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MUMASS;
            ROOT::Math::PtEtaPhiMVector a( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
        }

        //Bd
        double bVprob = TMath::Prob( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BdMassRangeTight[0] || svtMass->at(iB)>BdMassRangeTight[1] ) continue;
        if( bVprob < vtxProbMin ) continue;
        if(tB.Pt() < bsPtMin) continue;

        // int PV = GetPVPointing(iB, tB);
        // if(PV<0) continue;
        if(GetCt2D(tB, B0MASS, iB) < ctMin) continue;
        if(GetCt2D(tB, B0MASS, iB) / GetCt2DErr(tB, B0MASS, iB) < ctSigmaMin) continue;

        if( bVprob < best ) continue;
        index = iB;
        best = bVprob;
    }
    return index;
}

// ========================================================================================
int PhisUtil::GetBestJpsi()
{
    int index = -1;
    double bestChi2 = 1e9;
    for( int i=0; i<nSVertices; ++i ){
        if((svtType->at(i)!=PDEnumString::svtJPsi) ) continue;
        if( abs(svtMass->at(i)-JPSIMASS) > jpsiMassWin) continue;

        if( svtChi2->at(i)>bestChi2 ) continue;
        index = i;
        bestChi2 = svtChi2->at(i);
    }
    return index;
}

// ========================================================================================
double PhisUtil::GetInvMass(int i1, int i2, double mass1, double mass2)
{
    double px1 = trkPx->at(i1);
    double px2 = trkPx->at(i2);
    double py1 = trkPy->at(i1);
    double py2 = trkPy->at(i2);
    double pz1 = trkPz->at(i1);
    double pz2 = trkPz->at(i2);

    double E1 = sqrt( pow(mass1,2) + pow(px1,2)+pow(py1,2)+pow(pz1,2) );
    double E2 = sqrt( pow(mass2,2) + pow(px2,2)+pow(py2,2)+pow(pz2,2) );

    double m = sqrt( pow( (E1+E2),2) - ( pow(px1+px2,2) + pow(py1+py2,2) + pow(pz1+pz2,2) ) );

    return m;
}

// =====================================================================================
int PhisUtil::GetMixStatus( uint genindex )
{
    int Code = genId->at( genindex );

    const vector<int>& aD = allDaughters(genindex);
    if( aD.size()>0 && genId->at(aD[0]) == -Code ) return 2;

    const vector<int>& aM = allMothers(genindex); 
    if( aM.size()>0 && genId->at(aM[0]) == -Code ) return 1;

    return 0;
}

// ========================================================================================
double PhisUtil::GetMuoPFiso (int iMuon)
{
    double PFIso = muoSumCPpt->at(iMuon)/muoPt->at(iMuon);
    double betaCorr = muoSumNHet->at(iMuon)+muoSumPHet->at(iMuon)-0.5*(muoSumPUpt->at(iMuon));
    betaCorr/=muoPt->at(iMuon);
    if(betaCorr>0) PFIso+=betaCorr;

    return PFIso;
}


// =====================================================================================
double PhisUtil::GetJetCharge(int iJet, double kappa)
{
    double QJet = 0;
    double ptJet = 0;

    vector<int> list = pfCandFromJet( iJet );

    for(int it:list){
       double pt = pfcPt->at(it);
       double eta = pfcEta->at(it);

       if(pt<0.2) continue;
       if(abs(eta)>2.5) continue;

       QJet += pfcCharge->at(it) * pow(pt, kappa);
       ptJet += pow(pt, kappa);
    }

    QJet /= ptJet;

    return QJet; 
}

// =====================================================================================
double PhisUtil::GetListCharge(vector<int> *list, double kappa)
{
    double Q = 0;
    double pt = 0;
    for(int it:*list){
       Q += trkCharge->at(it) * pow(trkPt->at(it), kappa);
       pt += pow(trkPt->at(it), kappa);
    }
    return Q/pt; 
}

// =====================================================================================
double PhisUtil::GetJetProbb(int iJet)
{
    double probb = 0;
    double probbb = 0;
    double problepb = 0;
    for(int iTag=0; iTag<nTags; ++iTag){
        if(tagJet->at(iTag) != iJet) continue;
        if(tagType->at(iTag) == PDEnumString::pfDeepFlavourJetTags_probb){
            probb = tagProb->at(iTag);
            continue;
        }
        if(tagType->at(iTag) == PDEnumString::pfDeepFlavourJetTags_probbb){
            probbb = tagProb->at(iTag);
            continue;
        }
        if(tagType->at(iTag) == PDEnumString::pfDeepFlavourJetTags_problepb){
            problepb = tagProb->at(iTag);
            break;
        }
    }
    return probb + probbb + problepb;
}

// =====================================================================================
double PhisUtil::CountEventsWithFit(TH1 *hist, TString process)
{
    //ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 10000 );

    double mean = BSMASS;
    if(process=="BsJPsiPhi")   mean=BSMASS;
    if(process=="BuJPsiK")     mean=BUMASS;
    if(process=="BdJPsiKx")    mean=B0MASS;
    if(process=="BdKxMuMu")    mean=B0MASS;

    double sigma = 0.015;

    TString sgnDef = "[1]*TMath::Gaus(x, [0], [4], true)";
    sgnDef +=       "+[2]*TMath::Gaus(x, [0], [5], true)";
    sgnDef +=       "+[3]*TMath::Gaus(x, [0], [6], true)";
    TString bkgDef = "[7]+[8]*x";

    TString funcDef = sgnDef + "+" + bkgDef;

    TF1 *func = new TF1("func", funcDef, hist->GetBinLowEdge(1), hist->GetBinLowEdge(hist->GetNbinsX()));

    //SIGNAL
    double limit = hist->GetEntries()*hist->GetBinWidth(1);

    func->SetParameter(0, mean);
    func->SetParameter(1, limit/3);
    func->SetParameter(2, limit/3);
    func->SetParameter(3, limit/3);
    func->SetParameter(4, sigma);
    func->SetParameter(5, sigma);
    func->SetParameter(6, sigma);
    func->SetParLimits(0, mean-sigma, mean+sigma);
    func->SetParLimits(1, 0, limit);
    func->SetParLimits(2, 0, limit);
    func->SetParLimits(3, 0, limit);
    func->SetParLimits(4, sigma/2, sigma*2);
    func->SetParLimits(5, sigma/2, sigma*2);
    func->SetParLimits(6, sigma/2, sigma*2);

    //BKG
    func->SetParameter(7, 1);
    func->SetParameter(8, 100);

    hist->Fit("func","MRLQ");

    TF1 *fit = hist->GetFunction("MRLSQ");

    double nEvt = fit->GetParameter(1);
    nEvt += fit->GetParameter(2);
    nEvt += fit->GetParameter(3);

    nEvt/=hist->GetBinWidth(1);

    return nEvt;
}

// =====================================================================================
int PhisUtil::GetTightCandidate(TString process)
{
    if(process=="BsJPsiPhi") return GetBestBstrangeTight();
    if(process=="BuJPsiK")   return GetBestBupTight();
    return -1;
}

// =====================================================================================
int PhisUtil::GetCandidate(TString process)
{
    if(process=="BsJPsiPhi") return GetBestBstrange();
    if(process=="BuJPsiK")   return GetBestBup();
    return -1;
}

// =====================================================================================
double PhisUtil::dZ(int itk, int iPV)
{
    return PDAnalyzerUtil::dZ(itk, pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV));
}

// =====================================================================================
double PhisUtil::dXYjet(int itk, int iPV, int iJet)
{
    return abs(dXY( itk, pvtX->at(iPV), pvtY->at(iPV) ))*dSign( itk, iJet, pvtX->at(iPV), pvtY->at(iPV) );
}

// =================================================================================================
void PhisUtil::PrintMotherChain(int iGen, std::ostream& out)
{
    out<<genId->at(iGen)<<" << ";
    const vector<int>& vM = allMothers(iGen);
    uint nmot = vM.size();
    if(nmot>1) for(uint im=0; im<nmot; ++im) if(genId->at(vM[im])!=21) out<<genId->at(vM[im])<<" ";
    if(nmot==1) PrintMotherChain(vM[0]);
    return;
}

// =================================================================================================
void PhisUtil::PrintDaughterTree(int iGen, const string & pre)
{
    PrintDaughterTree(iGen, cout, pre);
}

// =================================================================================================
void PhisUtil::PrintDaughterTreePt(int iGen, const string & pre)
{
    PrintDaughterTree(iGen, cout, pre);
}

// =================================================================================================
void PhisUtil::PrintDaughterTree(int iGen, std::ostream& out, const std::string& indent, bool isLastSibling) {
  auto& vD = allDaughters(iGen);
  auto nDau = vD.size();
  
  std::string arrow = "|--->";
  if (nDau > 0) {
    arrow = "|-+->";
  }
  if (isLastSibling) {
    arrow[0] = '\\';
  }
  
  std::string newIndent = indent + isLastSibling?"  ":"| ";
  
  out << indent << arrow 
    << " pid: " << genId->at(iGen) 
    << " pt: " << genPt->at(iGen)
    << " eta: " << genEta->at(iGen) 
    << " phi: " << genPhi->at(iGen) 
    << std::endl;
  
  for (size_t iDau = 0; iDau < nDau; iDau++) {
    auto&& dau = vD[iDau];
    PrintDaughterTree(dau, out, newIndent, iDau + 1 == nDau);
  }
}

// =================================================================================================
bool PhisUtil::HasDaughter(int iGen)
{
    const vector<int>& vD = allDaughters(iGen);
    return vD.size()>0 ? true : false;
}
// double PhisUtil::GetCt3DErr(ROOT::Math::XYZVector vBs, int iSV, int iPV)
// {
//     using namespace ROOT::Math;
//     XYZPoint vSVT( svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV) );
//     XYZPoint vPV( pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV) );

//     auto vPointing = vSVT - vPV;

//     TMatrixD covSV(3,3);
//     double covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
//                       svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), 
//                       svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
//     covSV.SetMatrixArray(covSVArray);

//     TMatrixD covPV(3,3);
//     double covPVArray[]={pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
//                       pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV), 
//                       pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)};
//     covPV.SetMatrixArray(covPVArray);

//     TMatrixD covTot= covSV+covPV;

//     double distArray[]={double(vPointing.X()),double(vPointing.Y()),double(vPointing.Z())};
//     TVectorD diff(3,distArray);

//     if ( diff.Norm2Sqr()==0) return -1.; //if the secondary vertex is exactly the same as PV 

//     return BSMASS/t.P() * sqrt(covTot.Similarity(diff)) / sqrt(diff.Norm2Sqr()); 
// }

// =================================================================================================
void PhisUtil::SetJpsiMuCut(bool ctCut = true) // change only cuts different between the two HLT
{
    SetMuPtMin(muPtMin_jpsimu);
    SetKaonPtMin(kaonPtMin_jpsimu);
    SetCtSigmaMin(ctSigmaMin_jpsimu);
    if(ctCut) SetCtMin(ctMin_jpsimu);
    else      SetCtMin(-999.);
}

// =================================================================================================
void PhisUtil::SetJpsiTrkTrkCut(bool ctCut = true) // change only cuts different between the two HLT
{
    SetMuPtMin(muPtMin_jpsitrktrk);
    SetKaonPtMin(kaonPtMin_jpsitrktrk);
    SetCtSigmaMin(ctSigmaMin_jpsitrktrk);
    if(ctCut) SetCtMin(ctMin_jpsitrktrk);
    else      SetCtMin(-999.);
}

// =================================================================================================
int PhisUtil::HLTMatch(int trackIdx, bool isMu){
  
  int       hltIdx = -1;  
  double    minDist = 0.1; 

  // type 5 = track, 2 = muon
  for(int k=0; k<nHLTObjects; k++){

    if(hltObjType->at(k) != 2 && isMu) continue; // check only matches between trigger muons and muons
    if(hltObjType->at(k) != 5 && !isMu) continue; // check only matches between trigger tracks and kaons
    
    double dR = deltaR(trkEta->at(trackIdx), trkPhi->at(trackIdx), hltEta->at(k), hltPhi->at(k));

    if(dR < minDist){
      minDist = dR;
      hltIdx = k;
    }
  }
  return hltIdx;
}

// =================================================================================================
ROOT::Math::PxPyPzEVector PhisUtil::svtP4( int iSV, bool fromRefit ) {
  auto&& trkSsB = tracksFromSV(iSV);
  
  ROOT::Math::PxPyPzEVector sv(0,0,0,0);

  for (auto iSvTrk: trkSsB) {
    auto iSvTvp = vpIndex(iSvTrk, iSV);
    auto iSvTip = tvpTip->at(iSvTvp);
    
    ROOT::Math::PxPyPzEVector trk;
    if (fromRefit) trk = ROOT::Math::PtEtaPhiMVector( tvpPt->at(iSvTvp), tvpEta->at(iSvTvp), tvpPhi->at(iSvTvp), tipMass->at(iSvTip) );
    else trk = ROOT::Math::PtEtaPhiMVector( trkPt->at(iSvTrk), trkEta->at(iSvTrk), trkPhi->at(iSvTrk), tipMass->at(iSvTip) );
    
    sv += trk;
  }

  return sv;
}
