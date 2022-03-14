/////////////////////////////////////////////////////////////////////////////////
//      
//      Collection of utility functions for the Bs->J/PsiPhi CPV analysis
//
//      Author: Alberto Bragagnolo (alberto.bragagnolo@cern.ch)
//
/////////////////////////////////////////////////////////////////////////////////

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
 
    TLorentzVector pGen;
    pGen.SetPtEtaPhiM( (double) genPt->at(genIndex), (double) genEta->at(genIndex),
         (double) genPhi->at(genIndex), (double) genMass->at(genIndex) );

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
    // TLorentzVector tJPsi(0,0,0,0);           // No more requiremebt on Jpsi momentum

    for( uint i=0; i<tkJpsi.size(); ++i ){
        int j = tkJpsi[i];
        if(trkPt->at(j) < muPtMin) return false;
        if(abs(trkEta->at(j)) > muEtaMax) return false;
        // TLorentzVector a;
        // a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), MUMASS );
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

    for( int iB=0; iB<nSVertices; ++iB ){

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
        TLorentzVector tB(0,0,0,0);

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            double m = KMASS;
            if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MUMASS;
            TLorentzVector a;
            a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
        }

        double bVprob = TMath::Prob( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BsMassRangeTight[0] || svtMass->at(iB)>BsMassRangeTight[1] ) continue;
        if( bVprob < vtxProbMin ) continue;
        if(tB.Pt() < bsPtMin) continue;

        // int iPV = GetBestPV(iB, tB);
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

        TLorentzVector tB(0,0,0,0);
        double KaonPt = 0.;

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            double m = KMASS;
            if( j == tkJpsi[0] || j == tkJpsi[1] ){ m = MUMASS; }else{ KaonPt = trkPt->at(j); }
            TLorentzVector a;
            a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
       }

       //Bu
        double bVprob = TMath::Prob( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BuMassRangeTight[0] || svtMass->at(iB)>BuMassRangeTight[1] ) continue;
        if( bVprob < vtxProbMin ) continue;
        if(tB.Pt() < bsPtMin) continue;
        if(KaonPt < 1.6) continue;

        // int iPV = GetBestPV(iB, tB);
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

        TLorentzVector tB(0,0,0,0);

        for( uint i=0; i<tkSsB.size(); ++i ){
            int j = tkSsB[i];
            double m = KXMASS;
            if( j == tkJpsi[0] || j == tkJpsi[1] ) m = MUMASS;
            TLorentzVector a;
            a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
            tB += a;
        }

        //Bd
        double bVprob = TMath::Prob( svtChi2->at(iB), svtNDOF->at(iB) );
        if( svtMass->at(iB)<BdMassRangeTight[0] || svtMass->at(iB)>BdMassRangeTight[1] ) continue;
        if( bVprob < vtxProbMin ) continue;
        if(tB.Pt() < bsPtMin) continue;

        // int PV = GetBestPV(iB, tB);
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
int PhisUtil::TagMixStatus( uint genindex )
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
int PhisUtil::GetBestPV(int isvt, TLorentzVector t)   // TO BE CHECKED
{
    int ssPV = -1;
    double bestCos = -1.;
    double maxDz = 5.;

    TVector3 vB(t.Px(),t.Py(),t.Pz());
    TVector3 vSVT( svtX->at(isvt), svtY->at(isvt), svtZ->at(isvt) );

    for( int i=0; i<nPVertices; ++i ){
       if(abs(svtZ->at(isvt) - pvtZ->at( i )) > maxDz ) continue;

       TVector3 vPV(pvtX->at( i ), pvtY->at( i ), pvtZ->at( i ) );
       TVector3 vPointing;
       vPointing = vSVT - vPV;
       double cos = vPointing.Unit() * vB.Unit();

       if(cos > bestCos ){
         bestCos = cos;
         ssPV = i;
       }
    }
    return ssPV;
}

// =====================================================================================
TLorentzVector PhisUtil::GetTLorentzVecFromJpsiX(int iSvt)
{
    int iJPsi = (daugSVsFromVtx(iSvt, PDEnumString::vtrCascade)).at(0);
    vector<int> tkJpsi = tracksFromSV(iJPsi);
    vector<int> tkSsB = tracksFromSV(iSvt);
    TLorentzVector t(0,0,0,0);

    for( uint i=0; i<tkSsB.size(); ++i ){
        int j = tkSsB.at(i);
        double m = KMASS;
        if( (j==tkJpsi.at(0)) || (j==tkJpsi.at(1)) ) m = MUMASS;
        TLorentzVector a;
        a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
        t += a;
    }

    return t;
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
void PhisUtil::PrintMotherChain(int iGen)
{
    cout<<genId->at(iGen)<<" << ";
    const vector<int>& vM = allMothers(iGen);
    uint nmot = vM.size();
    if(nmot>1) for(uint im=0; im<nmot; ++im) if(genId->at(vM[im])!=21) cout<<genId->at(vM[im])<<" ";
    if(nmot==1) PrintMotherChain(vM[0]);
    return;
}

// =================================================================================================
void PhisUtil::PrintDaughterTree(int iGen, const string & pre)
{
    cout<<genId->at(iGen)<<endl;

    const vector<int>& vD = allDaughters(iGen);
    uint ndau = vD.size();
    if(ndau == 0) return;

    bool lastLevel = true;
    for(uint id =0; id<ndau; ++id){
        if ( HasDaughter( vD[id] ) ) {
            lastLevel = false;
            break;
        }
    }

    if( lastLevel ){
        cout<<pre<< "+-> ";
        for( uint id=0; id<ndau; ++id ) {
            int d = vD[id];
            cout<<genId->at(d)<<" ";
        }
        cout<<endl;
        return;
    }

    for( uint id=0; id<ndau; ++id ) {
        int d = vD[id];
        cout<<pre<< "+-> ";
        string prepre(pre);
        if ( id == ndau - 1 ) prepre += "    ";
        else prepre += "|   ";
        PrintDaughterTree( d, prepre );
    }
}

// =================================================================================================
void PhisUtil::PrintDaughterTreePt(int iGen, const string & pre)
{
    cout<<genId->at(iGen)<<" ["<<genPt->at(iGen)<<"]"<<endl;

    const vector<int>& vD = allDaughters(iGen);
    uint ndau = vD.size();
    if(ndau == 0) return;

    bool lastLevel = true;
    for(uint id =0; id<ndau; ++id){
        if ( HasDaughter( vD[id] ) ) {
            lastLevel = false;
            break;
        }
    }

    if( lastLevel ){
        cout<<pre<< "+-> ";
        for( uint id=0; id<ndau; ++id ) {
            int d = vD[id];
            cout<<genId->at(d)<<" ["<<genPt->at(iGen)<<"] ";
        }
        cout<<endl;
        return;
    }

  for( uint id=0; id<ndau; ++id ) {
    int d = vD[id];
    cout<<pre<< "+-> ";
    string prepre(pre);
    if ( id == ndau - 1 ) prepre += "    ";
    else prepre += "|   ";
    PrintDaughterTreePt( d, prepre );
  }
}

// =================================================================================================
bool PhisUtil::HasDaughter(int iGen)
{
    const vector<int>& vD = allDaughters(iGen);
    return vD.size()>0 ? true : false;
}

// =================================================================================================
double PhisUtil::GetCt2D(TLorentzVector t, double bMass, int iSV)
{
    TVector3 SVpos, PVpos;
    SVpos.SetXYZ(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));
    PVpos.SetXYZ(bsX, bsY, 0.);
    TVector3 Lxy = SVpos - PVpos;

    TVector3 ptB = t.Vect();
    ptB.SetZ(0.);

    return bMass*Lxy.Dot(ptB)/ptB.Mag2();
}

// =================================================================================================
double PhisUtil::GetCt2D(TLorentzVector t, double bMass, int iSV, int iPV)
{
    TVector3 SVpos, PVpos;
    SVpos.SetXYZ(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));
    PVpos.SetXYZ(pvtX->at(iPV),pvtY->at(iPV),0.);
    TVector3 Lxy = SVpos - PVpos;

    TVector3 ptB = t.Vect();
    ptB.SetZ(0.);
    
    return bMass*Lxy.Dot(ptB)/ptB.Mag2();
}

// =================================================================================================
double PhisUtil::GetCt3D(TLorentzVector t, double bMass, int iSV, int iPV)
{
    TVector3 SVpos, PVpos;
    SVpos.SetXYZ(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));
    PVpos.SetXYZ(pvtX->at(iPV),pvtY->at(iPV),pvtZ->at(iPV));
    TVector3 Lxyz = SVpos - PVpos;

    TVector3 pB = t.Vect();
    
    return bMass*Lxyz.Dot(pB)/pB.Mag2();
}

// =================================================================================================
double PhisUtil::GetCt2DErr(TLorentzVector t, double bMass, int iSV, int iPV)
{

    TVector3 SVpos, PVpos;
    SVpos.SetXYZ(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));
    PVpos.SetXYZ(pvtX->at(iPV),pvtY->at(iPV),0.);
    TVector3 Lxy = SVpos - PVpos;

    TVector3 ptB = t.Vect();
    ptB.SetZ(0.);

    double ct = bMass*Lxy.Dot(ptB)/ptB.Mag2();

    // Conversions for matrix multiplications
    double covSVArray[] = {
        svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
        svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV),
        svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)
    };

    double covPVArray[] = {
        pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
        pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV),
        pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)
    };

    TMatrixD covSV(3,3);
    TMatrixD covPV(3,3);
    covSV.SetMatrixArray(covSVArray);
    covPV.SetMatrixArray(covPVArray);

    TMatrixD covTot = covSV + covPV;

    double LxyArray[] = {Lxy.X(), Lxy.Y(), Lxy.Z()};
    TVectorD LxyVectorD(3, LxyArray);

    double ptBArray[] = {ptB.X(), ptB.Y(), ptB.Z()};
    TVectorD ptBVector(3, ptBArray);

    //Lxy contribution
    double ctErr_Lxy2 = pow(bMass/ptB.Mag2(),2)*covTot.Similarity(ptBVector);

    //Pt contribution
    vector<int> tkB  = tracksFromSV(iSV);
    int ntk = tkB.size();

    unordered_map<int,int> trkTotpp;

    for(int itpp=0; itpp<nTrkPer; itpp++)
        trkTotpp[tppTrk->at(itpp)] = itpp;

    TVector3 trkTVector3[ntk];
    int qArray[ntk];
    int trk_idx[ntk];

    for( int i=0; i<ntk; ++i ){
        int j = tkB.at(i);
        TVector3 tempTV3;
        tempTV3.SetXYZ(trkPx->at(j),trkPy->at(j),trkPz->at(j));

        trkTVector3[i]  = tempTV3;
        qArray[i]       = trkCharge->at(j);
        trk_idx[i]      = j;
    }

    double ctErr_pt2 = 0.;

    // Construct trkTotpp map
    for(int i=0; i<ntk; ++i){
        int itpp = trkTotpp.at(trk_idx[i]);

        double covHelixArray[] = {
            tppSQopQop->at(itpp), tppSQopLam->at(itpp), tppSQopPhi->at(itpp),
            tppSQopLam->at(itpp), tppSLamLam->at(itpp), tppSLamPhi->at(itpp),
            tppSQopPhi->at(itpp), tppSLamPhi->at(itpp), tppSPhiPhi->at(itpp)
        };
        TMatrixD covHelix(3,3);
        covHelix.SetMatrixArray(covHelixArray);

        // Helix parameters
        double Qop = qArray[i]/trkTVector3[i].Mag();
        double Lam = TMath::Pi()/2 - trkTVector3[i].Theta();

        double DQop = 1./(ptB.Mag2()*Qop)*(2*ct*ptB.Dot(trkTVector3[i]) - bMass*Lxy.Dot(trkTVector3[i]));
        double DLam = tan(Lam)/ptB.Mag2()*(2*ct*ptB.Dot(trkTVector3[i]) - bMass*Lxy.Dot(trkTVector3[i]));
        double DPhi = 1./ptB.Mag2()*(bMass*(Lxy.Y()*trkTVector3[i].X()-Lxy.X()*trkTVector3[i].Y())
                                  -2*ct*(ptB.Y()*trkTVector3[i].X()-ptB.X()*trkTVector3[i].Y()));

        double DHelixArray[] = {DQop, DLam, DPhi};
        TVectorD DHelixVectorD(3, DHelixArray);

        ctErr_pt2 += covHelix.Similarity(DHelixVectorD);
    }
    return sqrt(ctErr_Lxy2 + ctErr_pt2);
}


// =================================================================================================
double PhisUtil::GetCt2DErr(TLorentzVector t, double bMass, int iSV)
{

    TVector3 SVpos, PVpos;
    SVpos.SetXYZ(svtX->at(iSV),svtY->at(iSV),svtZ->at(iSV));
    PVpos.SetXYZ(bsX,bsY,0.);
    TVector3 Lxy = SVpos - PVpos;

    TVector3 ptB = t.Vect();
    ptB.SetZ(0.);

    double ct = bMass*Lxy.Dot(ptB)/ptB.Mag2();

    // Conversions for matrix multiplications
    double covSVArray[] = {
        svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
        svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV),
        svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)
    };

    double covPVArray[] = { // to be defined !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        0.,0.,0.,
        0.,0.,0.,
        0.,0.,0.
    };

    TMatrixD covSV(3,3);
    TMatrixD covPV(3,3);
    covSV.SetMatrixArray(covSVArray);
    covPV.SetMatrixArray(covPVArray);

    TMatrixD covTot = covSV + covPV;

    double LxyArray[] = {Lxy.X(), Lxy.Y(), Lxy.Z()};
    TVectorD LxyVectorD(3, LxyArray);

    double ptBArray[] = {ptB.X(), ptB.Y(), ptB.Z()};
    TVectorD ptBVector(3, ptBArray);

    //Lxy contribution
    double ctErr_Lxy2 = pow(bMass/ptB.Mag2(),2)*covTot.Similarity(ptBVector);

    //Pt contribution
    vector<int> tkB  = tracksFromSV(iSV);
    int ntk = tkB.size();

    unordered_map<int,int> trkTotpp;

    for(int itpp=0; itpp<nTrkPer; itpp++)
        trkTotpp[tppTrk->at(itpp)] = itpp;

    TVector3 trkTVector3[ntk];
    int qArray[ntk];
    int trk_idx[ntk];

    for( int i=0; i<ntk; ++i ){
        int j = tkB.at(i);
        TVector3 tempTV3;
        tempTV3.SetXYZ(trkPx->at(j),trkPy->at(j),trkPz->at(j));

        trkTVector3[i]  = tempTV3;
        qArray[i]       = trkCharge->at(j);
        trk_idx[i]      = j;
    }

    double ctErr_pt2 = 0.;

    // Construct trkTotpp map
    for(int i=0; i<ntk; ++i){
        int itpp = trkTotpp.at(trk_idx[i]);

        double covHelixArray[] = {
            tppSQopQop->at(itpp), tppSQopLam->at(itpp), tppSQopPhi->at(itpp),
            tppSQopLam->at(itpp), tppSLamLam->at(itpp), tppSLamPhi->at(itpp),
            tppSQopPhi->at(itpp), tppSLamPhi->at(itpp), tppSPhiPhi->at(itpp)
        };
        TMatrixD covHelix(3,3);
        covHelix.SetMatrixArray(covHelixArray);

        // Helix parameters
        double Qop = qArray[i]/trkTVector3[i].Mag();
        double Lam = TMath::Pi()/2 - trkTVector3[i].Theta();

        double DQop = 1./(ptB.Mag2()*Qop)*(2*ct*ptB.Dot(trkTVector3[i]) - bMass*Lxy.Dot(trkTVector3[i]));
        double DLam = tan(Lam)/ptB.Mag2()*(2*ct*ptB.Dot(trkTVector3[i]) - bMass*Lxy.Dot(trkTVector3[i]));
        double DPhi = 1./ptB.Mag2()*(bMass*(Lxy.Y()*trkTVector3[i].X()-Lxy.X()*trkTVector3[i].Y())
                                  -2*ct*(ptB.Y()*trkTVector3[i].X()-ptB.X()*trkTVector3[i].Y()));

        double DHelixArray[] = {DQop, DLam, DPhi};
        TVectorD DHelixVectorD(3, DHelixArray);

        ctErr_pt2 += covHelix.Similarity(DHelixVectorD);
    }
    
    return sqrt(ctErr_Lxy2 + ctErr_pt2);
}


// =================================================================================================
// double PhisUtil::GetCt3DErr(TLorentzVector t, int iSV, int iPV)
// {
//     TVector3 vSVT( svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV) );
//     TVector3 vPV( pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV) );

//     TVector3 vPointing = vSVT - vPV;
//     TVector3 vBs = t.Vect();

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
void PhisUtil::SetJpsiMuCut(bool ctCut = true) // change only different cuts
{
    SetBsPtMin(10);
    SetMuPtMin(3.5);
    SetKaonPtMin(1.2);
    if(ctCut) SetCtMin(0.007);
    if(ctCut) SetCtSigmaMin(-999.);       // no cut
}

// =================================================================================================
void PhisUtil::SetJpsiTrkTrkCut(bool ctCut = true) // change only different cuts
{
    SetBsPtMin(9);
    SetMuPtMin(4.0);
    SetKaonPtMin(0.9);
    if(ctCut) SetCtMin(-999.);            // no cut
    if(ctCut) SetCtSigmaMin(3);
}
