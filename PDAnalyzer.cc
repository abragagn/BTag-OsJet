#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>
#include <algorithm>

#include "PDAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TFile.h"

// additional features
#include "PDSecondNtupleWriter.h"   // second ntuple
//#include "DataSetFilter.cc"   // dataset filter
#include "PDMuonVar.cc"
#include "PDSoftMuonMvaEstimator.cc"
#include "AlbertoUtil.cc"
#include "OSMuonMvaTag.cc"

using namespace std;

/*
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list histHevjin.root -v outputFile ntuHevjin.root -v histoMode RECREATE -v use_gen t -n 10000
*/

PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // user parameters are set as names associated to a string, 
    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useHLT", "true" );

    //Jet Parameters
    setUserParameter( "CutDeepCSV", "0.4941" );  //loose = 0.1522, medium = 0.4941, tight = 0.8001
    setUserParameter( "kappa", "1.0" );
    setUserParameter( "QCut", "0.0" );
    setUserParameter( "minPtJet", "10" );
    setUserParameter( "jetDrCut", "0.5" );
    setUserParameter( "jetDzCut", "0.2" );
    setUserParameter( "jetChargeDzCut", "0.2" );

    setUserParameter( "ptCut", "40.0" ); //needed for paolo's code for unknow reasons

}


PDAnalyzer::~PDAnalyzer() {
}



void PDAnalyzer::beginJob() {

    PDAnalyzerUtil::beginJob();

    // user parameters are retrieved as strings by using their names;
    // numeric parameters ( int, float or whatever ) can be directly set
    // by passing the corresponding variable,
    // e.g. getUserParameter( "name", x )

    getUserParameter( "verbose", verbose );

    getUserParameter( "useHLT", useHLT );
    getUserParameter( "process", process );
 
    getUserParameter( "minPtJet", minPtJet );
    getUserParameter( "kappa", kappa );
    getUserParameter( "CutDeepCSV", CutDeepCSV );
    getUserParameter( "QCut", QCut );
    getUserParameter( "minPtJet", minPtJet );
    getUserParameter( "jetDrCut", jetDrCut );
    getUserParameter( "jetDzCut", jetDzCut );
    getUserParameter( "jetChargeDzCut", jetChargeDzCut );

    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

// additional features

    tWriter = new PDSecondNtupleWriter; // second ntuple
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    if(process=="BsJPsiPhi") SetBsMassRange(5.20, 5.50);
    if(process=="BuJPsiK") SetBuMassRange(5.1, 5.50);

    setOsMuonCuts(0.8910, 0.8925, 1. ); //set wp barrel, wp endcap for muonID and Dz cut
    inizializeMuonMvaReader( "BDTMuonID2017woIPwIso" ); //initialize mva muon id reader
    inizializeOSMuonMvaTagReader( "DNNOsMuonHLTJpsiMu_test241" ); //inizialize mva os muon reader
    inizializeOSMuonMvaMistagMethods();   //read os muon method for per-event-mistag 

    return;

}


void PDAnalyzer::book() {

    // putting "autoSavedObject" in front of the histo creation 
    // it's automatically marked for saving on file; the option 
    // is uneffective when not using the full utility

    return;

}


void PDAnalyzer::reset() {
// automatic reset
    autoReset();
    return;
}


bool PDAnalyzer::analyze( int entry, int event_file, int event_tot ) {

    if ( verbose ) {
        cout << " +++++++++++++++++++++++++++ " << endl;
        cout << "entry: "
             << entry << " " << event_file << " " << event_tot << endl;
        cout << "run: " << runNumber << " , "
             << "evt: " << eventNumber << endl;
    }
    else {

        if ( (!(event_tot%10) && event_tot<100 ) || 
     (!(event_tot %100) && event_tot<1000 ) || 
     (!(event_tot %1000) && event_tot<10000 ) || 
     (!(event_tot %10000) && event_tot<100000 ) || 
     (!(event_tot %100000) && event_tot<1000000 ) || 
     (!(event_tot %1000000) && event_tot<10000000 ) )
            cout << " == at event " << event_file << " " << event_tot << endl;
    }
// additional features
    computeMuonVar();   //compute muon variable for soft id
    inizializeTagVariables(); //initialiaze some variable for tagging
    tWriter->Reset();
    convSpheCart(jetPt, jetEta, jetPhi, jetPx, jetPy, jetPz);
    convSpheCart(muoPt, muoEta, muoPhi, muoPx, muoPy, muoPz);
    convSpheCart(trkPt, trkEta, trkPhi, trkPx, trkPy, trkPz);
    convSpheCart(pfcPt, pfcEta, pfcPhi, pfcPx, pfcPy, pfcPz);

    if( !((process=="BsJPsiPhi")||(process=="BuJPsiK")) ) {
        cout<<"!$!#$@$% PROCESS NAME WRONG"<<endl;
        return false;
    }

//------------------------------------------------HLT---------------------------------------

    bool jpsimu = false;
    bool jpsitktk = false;
    bool jpsitk = false;

    if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v)||hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v)) jpsimu = true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v)) jpsitktk =  true;
    if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v)) jpsitk = true;

    if( jpsimu ) return false;
    SetJpsiTrkTrkCut();

    if(useHLT && process=="BsJPsiPhi" && !jpsitktk) return false;
    if(useHLT && process=="BuJPsiK" && !jpsitk) return false;

//------------------------------------------------SEARCH FOR SS---------------------------------------

    int ssbSVT = GetCandidate(process);
    if(ssbSVT<0) return false;

    bool isTight = false;
    int ssbSVTtight = GetTightCandidate(process);
    if(ssbSVTtight>=0){
        isTight = true;
        ssbSVT = ssbSVTtight;
    }

    int iJPsi = (subVtxFromSV(ssbSVT)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi);
    vector <int> tkSsB = tracksFromSV(ssbSVT);

    TLorentzVector tB = GetTLorentzVecFromJpsiX(ssbSVT);

    //generation information

    vector <int> ListLongLivedB;
    vector <int> ListB;
    int genBindex = -1;
    int ssbLund = 0;
    int tagMix = -1;
    float evtWeight = 1;

    if(use_gen){
        for( unsigned int i=0 ; i<genId->size() ; ++i ){
            if(TagMixStatus( i ) == 2) continue;
            if( IsB(i) ) ListB.push_back(i);
            unsigned int Code = abs(genId->at(i));
            if( Code == 511 || Code == 521 || Code == 531 || Code == 541 || Code == 5122 ) ListLongLivedB.push_back(i);
        }

        genBindex = GetClosestGenLongLivedB( tB.Eta(), tB.Phi(), tB.Pt(), &ListLongLivedB);
        if(genBindex<0) return false;

        ssbLund = genId->at(genBindex);
        if((process=="BsJPsiPhi") && (abs(ssbLund)!=531)) return false;
        if((process=="BuJPsiK") && (abs(ssbLund)!=521)) return false;

        tagMix = TagMixStatus( genBindex );
        if(tagMix == 2) return false;
        if(tagMix == 1) ssbLund*=-1;

        for(auto it:ListLongLivedB){
            if(it == genBindex) continue;
            if(abs(genId->at(it)) == abs(ssbLund)) evtWeight = 2;
        }


    }else{
        if(process=="BsJPsiPhi") ssbLund = ((double)rand() / (RAND_MAX)) < 0.5 ? +531 : -531; //this should not be used
        if(process=="BuJPsiK"){
            for( auto it:tkSsB ){
                if( it == tkJpsi[0] || it == tkJpsi[1] ) continue;
                ssbLund = trkCharge->at(it) > 0 ? +521 : -521;
            }
        }
    }

    int ssbPVT = GetBestPV(ssbSVT, tB);
    if(ssbPVT < 0) return false;
    setSsForTag(ssbSVT, ssbPVT);

    if(getOsMuon()>=0) return false; 
    if(nElectrons>0) return false; 

    //FILLING SS
    (tWriter->ssbMass) = svtMass->at(ssbSVT);
    (tWriter->ssbIsTight) = isTight;

    (tWriter->ssbSVT) = ssbSVT;
    (tWriter->ssbPVT) = ssbPVT;
    (tWriter->ssbPVTx) = pvtX->at(ssbPVT);
    (tWriter->ssbPVTy) = pvtY->at(ssbPVT);
    (tWriter->ssbPVTz) = pvtZ->at(ssbPVT);
    
    (tWriter->ssbLund) = ssbLund;

    (tWriter->hltJpsiMu) = jpsimu;
    (tWriter->hltJpsiTrkTrk) = jpsitktk;
    (tWriter->hltJpsiTrk) = jpsitk;
    
    (tWriter->evtWeight) = evtWeight;

//-----------------------------------------JET-----------------------------------------

    int bestJet = -1;
    float bestJetTag = CutDeepCSV;

    //SELECTION
    for (int iJet = 0; iJet<nJets; ++iJet){

        vector <int> jet_pfcs = pfCandFromJet( iJet );
        vector <int> jet_tks = tracksFromJet( iJet );

        //SELECTION
        if(goodJet(iJet)!=true) continue;
        if(abs(jetEta->at(iJet))>2.5) continue;
        if(jetPt->at(iJet)<minPtJet) continue;

        float bTag = GetJetProbb(iJet);
        float cutbTag = CutDeepCSV;
        if(bTag < cutbTag) continue;

        float jetDrB = deltaR(jetEta->at( iJet ), jetPhi->at( iJet ), tB.Eta(), tB.Phi());
        if( jetDrB < jetDrCut ) continue;

        bool skip = false; 
        for(auto it:jet_tks){
            if(std::find(tkSsB.begin(), tkSsB.end(), it) != tkSsB.end()){
                skip = true;
                break;
            }
        }

        if(skip) continue;

        int nTrkNear = 0;
        for(auto it:jet_tks){
            if(fabs(dZ(it, ssbPVT))>=jetDzCut) continue;
            if( !(( trkQuality->at( it ) >> 2 ) & 1) ) continue;
            if( fabs(trkEta->at(it))>2.5 ) continue;
            if( trkPt->at(it)<0.7 ) continue;
            nTrkNear++;
        }
        if(nTrkNear < 2) continue;

        if(bTag>bestJetTag){
            bestJetTag = bTag;
            bestJet = iJet;
        }
    }

    //TAG
    if(bestJet < 0) return true;

    //indices
    int iJet = bestJet;
    vector <int> jet_tks = tracksFromJet( iJet );

    vector<int> selectedJetTracks;
    for(auto it:jet_tks){
        if(fabs(dZ(it, ssbPVT))>=jetChargeDzCut) continue;
        if( !(( trkQuality->at( it ) >> 2 ) & 1) ) continue;
        if( fabs(trkEta->at(it))>2.5 ) continue;
        if( trkPt->at(it)<0.7 ) continue;
        selectedJetTracks.push_back(it);
    }

    if(selectedJetTracks.size()==0) return true;


    //GENINFO
    int jetAncestor = GetJetAncestor( iJet, &ListB );

    //TAGGING VARIABLES
    //Separation
    float jetDzB=0;
    for(auto it:selectedJetTracks)
        jetDzB += fabs(dZ(it, ssbPVT));

    jetDzB = jetDzB/=selectedJetTracks.size();
    float jetDrB = deltaR(jetEta->at( iJet ), jetPhi->at( iJet ), tB.Eta(), tB.Phi());

    //Jet Charge
    float jet_charge = GetListCharge(&selectedJetTracks, kappa);

    (tWriter->jetPt) = jetPt->at(iJet);
    (tWriter->jetEta) = jetEta->at(iJet);
    (tWriter->jetPhi) = jetPhi->at(iJet);
    (tWriter->jetCharge) = jet_charge;
    (tWriter->jetDeepCSV) = GetJetProbb(iJet);
    (tWriter->jetDrB) = jetDrB;
    (tWriter->jetDzB) = jetDzB;
    (tWriter->jetNDau) = jetNDau->at(iJet);
    (tWriter->jetNHF) = jetNHF->at(iJet);
    (tWriter->jetNEF) = jetNEF->at(iJet);
    (tWriter->jetCHF) = jetCHF->at(iJet);
    (tWriter->jetCEF) = jetCEF->at(iJet);
    (tWriter->jetNCH) = jetNCH->at(iJet);

    //------------------------------------------------TAG------------------------------------------------

    (tWriter->jetHasAncestor) = jetAncestor;

    //------------------------------------------------TRACKS------------------------------------------------

    for(auto it:jet_tks){
        (tWriter->trkPt)->push_back(trkPt->at(it));
        (tWriter->trkEta)->push_back(trkEta->at(it));
        (tWriter->trkPhi)->push_back(trkPhi->at(it));
        (tWriter->trkCharge)->push_back(trkCharge->at(it));
        (tWriter->trkIsHighPurity)->push_back( ( trkQuality->at( it ) >> 2 ) & 1 );
        (tWriter->trkDxy)->push_back( dXYjet(it, ssbPVT, iJet) );
        (tWriter->trkDz)->push_back(dZ(it, ssbPVT));        
        (tWriter->trkIsInJet)->push_back( 1 );
    }

    for (int iTrk = 0; iTrk<nTracks; ++iTrk){

        if( deltaR(jetEta->at(iJet), jetPhi->at(iJet), trkEta->at(iTrk), trkPhi->at(iTrk)) > 0.5 ) continue;
        if( fabs(dZ(iTrk, ssbPVT)) > 0.1 ) continue;
        if(std::find(jet_tks.begin(), jet_tks.end(), iTrk) != jet_tks.end()) continue;

        (tWriter->trkPt)->push_back(trkPt->at(iTrk));
        (tWriter->trkEta)->push_back(trkEta->at(iTrk));
        (tWriter->trkPhi)->push_back(trkPhi->at(iTrk));
        (tWriter->trkCharge)->push_back(trkCharge->at(iTrk));
        (tWriter->trkIsHighPurity)->push_back( ( trkQuality->at( iTrk ) >> 2 ) & 1 );
        (tWriter->trkDxy)->push_back( dXYjet(iTrk, ssbPVT, iJet) );
        (tWriter->trkDz)->push_back(dZ(iTrk, ssbPVT));
        (tWriter->trkIsInJet)->push_back( 0 );        
    }

    (tWriter->evtNumber) = event_tot;
    tWriter->fill();

    return true;

}


void PDAnalyzer::endJob() {


// additional features

    tWriter->close();   // second ntuple

    return;
}


void PDAnalyzer::save() {
#   if UTIL_USE == FULL
    // explicit saving not necessary for "autoSavedObjects"
    autoSave();
#elif UTIL_USE == BARE
    // explicit save histos when not using the full utility

#endif

    return;
}


// ======MY FUNCTIONS===============================================================================
// =====================================================================================
int PDAnalyzer::GetJetAncestor( unsigned int iJet, vector<int> *GenList )
{
    double drb = 0.4;
    double dpb = minPtJet; 
    int best = -1;

    for(auto it:*GenList){

        float dr = deltaR(jetEta->at(iJet), jetPhi->at(iJet), genEta->at(it), genPhi->at(it));
        float pt = genPt->at(it);

        if( dr > drb ) continue;
        if( pt < dpb) continue;

        best = it;
        dpb = pt;

    }

    return best;
}
// =====================================================================================
int PDAnalyzer::GetClosesJet( float pt, float eta, float phi )
{
    double drb = 0.5;
    int best = -1;
    for( int i=0; i<nJets; ++i ){
       float dr = deltaR(eta, phi, jetEta->at(i), jetPhi->at(i));
       if( dr > drb ) continue;
       best = (int) i;
       drb = dr;
    }
    return best;
}
