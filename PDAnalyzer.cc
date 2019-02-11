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

using namespace std;

//pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2016Lists/BsToJpsiPhi_BMuonFilter_AOD_DCAP.list hist.root -v histoMode RECREATE -v use_gen t -n 10000

PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // user parameters are set as names associated to a string, 
    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useHLT", "false" );

    //Jet Parameters
    setUserParameter( "useJet", "t" ); 
    setUserParameter( "CutCSV", "0.8484" );
    setUserParameter( "CutDF", "0.5" );
    setUserParameter( "kappa", "0.8" );
    setUserParameter( "QCut", "0.0" );
    setUserParameter( "minPtJet", "16" );
    setUserParameter( "jetSeparationCut", "0.4" );
    setUserParameter( "jetDzCut", "10" );

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
    getUserParameter( "CutCSV", CutCSV );
    getUserParameter( "CutDF", CutDF );
    getUserParameter( "QCut", QCut );
    getUserParameter( "minPtJet", minPtJet );
    getUserParameter( "jetSeparationCut", jetSeparationCut );
    getUserParameter( "jetDzCut", jetDzCut );

    getUserParameter( "ptCut", ptCut ); //needed for paolo's code for unknow reasons

// additional features

    tWriter = new PDSecondNtupleWriter; // second ntuple
    tWriter->open( getUserParameter("outputFile"), "RECREATE" ); // second ntuple

    if(process=="BsJPsiPhi") SetBsMassRange(5.20, 5.50);
    if(process=="BuJPsiK") SetBuMassRange(5.1, 5.50);

    return;

}


void PDAnalyzer::book() {

    // putting "autoSavedObject" in front of the histo creation 
    // it's automatically marked for saving on file; the option 
    // is uneffective when not using the full utility


    float min = 5.0;
    float max = 6.0;
    float nbin = 500;


    autoSavedObject =
    hmass_ssB       = new TH1D( "hmass_ssB", "hmass_ssB", nbin, min, max );

    autoSavedObject =
    hmass_ssB_os        = new TH1D( "hmass_ssB_os", "hmass_ssB_os", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osWT  = new TH1D( "hmass_ssB_osWT", "hmass_ssB_osWT", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osRT  = new TH1D( "hmass_ssB_osRT", "hmass_ssB_osRT", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osJetwAnc = new TH1D( "hmass_ssB_osJetwAnc", "hmass_ssB_osJetwAnc", nbin, min, max );
    
    autoSavedObject =
    hmass_ssB_osJetwoAnc    = new TH1D( "hmass_ssB_osJetwoAnc", "hmass_ssB_osJetwoAnc", nbin, min, max );

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
    tWriter->Reset();
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

    if( !jpsitktk ) return false;
    SetJpsiTrkTrkCut();

    if(useHLT && process=="BsJPsiPhi" && !(jpsimu || jpsitktk)) return false;
    if(useHLT && process=="BuJPsiK" && !(jpsimu || jpsitk)) return false;

//------------------------------------------------SEARCH FOR SS---------------------------------------

    int iSsB = GetCandidate(process);
    if(iSsB<0) return false;

    bool isTight = false;
    int iSsBtight = GetTightCandidate(process);
    if(iSsBtight>=0){
        isTight = true;
        iSsB = iSsBtight;
    }

    int iJPsi = (subVtxFromSV(iSsB)).at(0);
    vector <int> tkJpsi = tracksFromSV(iJPsi);
    vector <int> tkSsB = tracksFromSV(iSsB);

    TLorentzVector tB = GetTLorentzVecFromJpsiX(iSsB);

    //generation information

    vector <int> ListLongLivedB;
    vector <int> ListB;
    int genBindex = -1;
    int ssBLund = 0;
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

        ssBLund = genId->at(genBindex);
        if((process=="BsJPsiPhi") && (abs(ssBLund)!=531)) return false;
        if((process=="BuJPsiK") && (abs(ssBLund)!=521)) return false;

        tagMix = TagMixStatus( genBindex );
        if(tagMix == 2) return false;
        if(tagMix == 1) ssBLund*=-1;

        for(auto it:ListLongLivedB){
            if(it == genBindex) continue;
            if(abs(genId->at(it)) == abs(ssBLund)) evtWeight = 2;
        }


    }else{
        if(process=="BsJPsiPhi") ssBLund = ((double)rand() / (RAND_MAX)) < 0.5 ? +531 : -531; //this should not be used
        if(process=="BuJPsiK"){
            for( auto it:tkSsB ){
                if( it == tkJpsi[0] || it == tkJpsi[1] ) continue;
                ssBLund = trkCharge->at(it) > 0 ? +521 : -521;
            }
        }
    }

    int iSsPV = GetBestPV(iSsB, tB);
    if(iSsPV < 0) return false;

    //FILLING SS
    (tWriter->ssbPt) = tB.Pt();
    (tWriter->ssbEta) = tB.Eta();
    (tWriter->ssbPhi) = tB.Phi();
    (tWriter->ssbMass) = svtMass->at(iSsB);
    (tWriter->ssbIsTight) = isTight;

    (tWriter->ssbLxy) = GetCt2D(tB, iSsB) / (MassBs/tB.Pt());
    (tWriter->ssbCt2D) = GetCt2D(tB, iSsB);
    (tWriter->ssbCt2DErr) = GetCt2DErr(tB, iSsB, iSsPV);
    (tWriter->ssbCt2DSigmaUnit) = GetCt2D(tB, iSsB, iSsPV)/GetCt2DErr(tB, iSsB, iSsPV);
    (tWriter->ssbCt3D) = GetCt3D(tB, iSsB, iSsPV);
    (tWriter->ssbCt3DErr) = GetCt3DErr(tB, iSsB, iSsPV);
    (tWriter->ssbCt3DSigmaUnit) = GetCt3D(tB, iSsB, iSsPV)/GetCt3DErr(tB, iSsB, iSsPV);

    (tWriter->ssbSVT) = iSsB;
    (tWriter->ssbPVT) = iSsPV;
    
    (tWriter->ssbLund) = ssBLund;

    (tWriter->hltJpsiMu) = jpsimu;
    (tWriter->hltJpsiTrkTrk) = jpsitktk;
    (tWriter->hltJpsiTrk) = jpsitk;
    
    (tWriter->evtWeight) = evtWeight;
    (tWriter->evtNb) = ListLongLivedB.size();

    hmass_ssB->Fill(svtMass->at(iSsB), evtWeight);

//-----------------------------------------OPPOSITE SIDE-----------------------------------------

    int bestJet = -1;
    float bestJetPt = minPtJet; 

    //SELECTION
    for (int iJet = 0; iJet<nJets; ++iJet){

        vector <int> jet_pfcs = pfCandFromJet( iJet );
        vector <int> jet_tks = tracksFromJet( iJet );

        //SELECTION
        if(goodJet(iJet)!=true) continue;
        if(jet_pfcs.size()<2) continue;
        if(abs(jetEta->at(iJet))>2.4) continue;
        if(jetPt->at(iJet)<minPtJet) continue;
        if(GetJetProbb(iJet)>=0){
            if(GetJetProbb(iJet)<CutDF) continue;
        }else{
            if(jetCSV->at(iJet)<CutCSV) continue;
        }

        float jetDrB = deltaR(jetEta->at( iJet ), jetPhi->at( iJet ), tB.Eta(), tB.Phi());
        if( jetDrB < jetSeparationCut ) continue;

        bool skip = false; 
        for(auto it:jet_tks)
            if(std::find(tkSsB.begin(), tkSsB.end(), it) != tkSsB.end()) 
                skip = true;

        if(skip) continue;

        float jetDz=0;
        for(auto it:jet_tks) jetDz += dZ(it, iSsPV);
        jetDz/=jet_tks.size();
        if(abs(jetDz)>10) continue;

        if(jetPt->at(iJet)>bestJetPt){
            bestJetPt = jetPt->at(iJet);
            bestJet = iJet;
        }

    }

    //TAG
    if(bestJet < 0){
        (tWriter->osJet) = 0;
        (tWriter->evtNumber)= event_tot;
        tWriter->fill();
        return true;
    }

    (tWriter->osJet) = 1;

    hmass_ssB_os->Fill(svtMass->at(iSsB));

    //indices
    int iJet = bestJet;
    vector <int> jet_pfcs = pfCandFromJet( iJet );
    vector <int> jet_tks = tracksFromJet( iJet );

    //GENINFO
    int jetAncestor = GetJetAncestor( iJet, &ListB );

    //TAGGING VARIABLES
    //Separation
    float jetDz=0;
    for(auto it:jet_tks) jetDz+= dZ(it, iSsPV);

    float jetDzB = jetDz/jet_tks.size();

    float jetDrB = deltaR(jetEta->at( iJet ), jetPhi->at( iJet ), tB.Eta(), tB.Phi());

    //Jet Charge
    float jet_charge = GetJetCharge(iJet, kappa);

    //------------------------------------------------FILLING------------------------------------------------
    (tWriter->jetPt)= jetPt->at(iJet);
    (tWriter->jetEta)= jetEta->at(iJet);
    (tWriter->jetPhi)= jetPhi->at(iJet);
    (tWriter->jetCharge)= jet_charge;
    (tWriter->jetCSV)= jetCSV->at(iJet);
    (tWriter->jetDFprobb)= GetJetProbb(iJet);
    (tWriter->jetSize)= jet_tks.size();
    (tWriter->jetDrB)= jetDrB;
    (tWriter->jetDzB)= jetDzB ;

    //------------------------------------------------TAG------------------------------------------------
    
    if(jetAncestor<0){
        hmass_ssB_osJetwoAnc->Fill(svtMass->at(iSsB));
        (tWriter->jetHasAncestor)= 0;

    }else{
        hmass_ssB_osJetwAnc->Fill(svtMass->at(iSsB));
        (tWriter->jetHasAncestor)= 1;
    }

    int isB = 0;
    (tWriter->osJetTag) = -2;

    if(jet_charge < -QCut ) isB = +1;
    if(jet_charge > +QCut ) isB = -1;

    if( TMath::Sign(1, ssBLund) == isB ){ 
        hmass_ssB_osRT->Fill(svtMass->at(iSsB));
        (tWriter->osJetTag)= 1;
    }

    if( TMath::Sign(1, ssBLund) == -1*isB ){
        hmass_ssB_osWT->Fill(svtMass->at(iSsB));
        (tWriter->osJetTag)= 0;
    }

    (tWriter->evtNumber)= event_tot;
    tWriter->fill();

    return true;

}


void PDAnalyzer::endJob() {


// additional features

    tWriter->close();   // second ntuple

//printing results

    float eff   = static_cast<float>(hmass_ssB_os->GetEntries()) / static_cast<float>(hmass_ssB->GetEntries());
    float w     = static_cast<float>(hmass_ssB_osWT->GetEntries()) / static_cast<float>(hmass_ssB_os->GetEntries());
    float power = eff*pow(1-2*w, 2);

    cout<<"wAnc woAnc"<<endl;
    cout<<hmass_ssB_osJetwAnc->GetEntries()<<" "<<hmass_ssB_osJetwoAnc->GetEntries()<<endl<<endl;

    cout<<"#Bp eff% w% P%"<<endl;
    cout<<hmass_ssB->GetEntries()<<" "<<eff*100<<" "<<w*100<<" "<<power*100<<endl;

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


// to plot some histogram immediately after the ntuple loop
// "uncomment" the following lines
/*
void PDAnalyzer::plot() {
    TCanvas* can = new TCanvas( "muoPt", "muoPt", 800, 600 );
    can->cd();
    can->Divide( 1, 2 );
    can->cd( 1 );
    hptmumax->Draw();
    hptmu2nd->Draw();
    return;
}
*/


// ======MY FUNCTIONS===============================================================================
// =====================================================================================
int PDAnalyzer::GetJetAncestor( unsigned int iJet, vector<int> *GenList ){

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
