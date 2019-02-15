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

/*
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list hist.root -v outputFile ntu.root -v histoMode RECREATE -v use_gen t -v useHLT t -n 10000
*/

PDAnalyzer::PDAnalyzer() {

    std::cout << "new PDAnalyzer" << std::endl;

    // user parameters are set as names associated to a string, 
    // default values can be set in the analyzer class contructor

    setUserParameter( "verbose", "f" );

    setUserParameter( "process", "BsJPsiPhi" );
    setUserParameter( "useHLT", "false" );

    setUserParameter( "saveNotTaggedEvets", "true" );

    //Jet Parameters
    setUserParameter( "CutCSV", "0.8838" ); //loose = 0.5803, medium = 0.8838, tight = 0.9693
    setUserParameter( "CutDeepCSV", "0.4941" );  //loose = 0.1522, medium = 0.4941, tight = 0.8001
    setUserParameter( "kappa", "1.0" );
    setUserParameter( "QCut", "0.0" );
    setUserParameter( "minPtJet", "14" );
    setUserParameter( "jetSeparationCut", "0.4" );
    setUserParameter( "jetDzCut", "5" );

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

    getUserParameter( "saveNotTaggedEvets", saveNotTaggedEvets );
 
    getUserParameter( "minPtJet", minPtJet );
    getUserParameter( "kappa", kappa );
    getUserParameter( "CutCSV", CutCSV );
    getUserParameter( "CutDeepCSV", CutDeepCSV );
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
    hmass_ssB       = new TH1F( "hmass_ssB", "hmass_ssB", nbin, min, max );

    autoSavedObject =
    hmass_ssB_os        = new TH1F( "hmass_ssB_os", "hmass_ssB_os", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osWT  = new TH1F( "hmass_ssB_osWT", "hmass_ssB_osWT", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osRT  = new TH1F( "hmass_ssB_osRT", "hmass_ssB_osRT", nbin, min, max );

    autoSavedObject =
    hmass_ssB_osJetwAnc = new TH1F( "hmass_ssB_osJetwAnc", "hmass_ssB_osJetwAnc", nbin, min, max );
    
    autoSavedObject =
    hmass_ssB_osJetwoAnc    = new TH1F( "hmass_ssB_osJetwoAnc", "hmass_ssB_osJetwoAnc", nbin, min, max );

    autoSavedObject =
    hTest    = new TH1F( "hTest", "hTest", 500, -2, 2 );

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

    if( jpsimu ) return false;
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

    vector <int> trk_pv = tracksFromPV( iSsPV );
    vector <int> trk_pv0 = tracksFromPV( 0 );

    vector<int> pvTrkSVT;
    for(auto it:tkSsB) pvTrkSVT.push_back(trkPVtx->at(it));
    int max = 0;
    int most_common = -1;
    map<int,int> m;
    for (auto vi:pvTrkSVT) {
      m[vi]++;
      if (m[vi] > max) {
        max = m[vi]; 
        most_common = vi;
      }
    }

    vector <int> trk_pvcommon = tracksFromPV( most_common );

/////////// DEBUG VERTICES
/*
    int selFlag = PDEnumString::general | PDEnumString::packed;
    vector<int> testList;

    cout<<endl<<endl<<"SVT "<<most_common<<endl<<endl;
    for(auto it:tkSsB){
        if(!(( trkQuality->at( it ) >> 2 ) & 1)) continue;
        cout<<it<<" - ";
        int gen = GetClosestGen(trkEta->at(it), trkPhi->at(it), trkPt->at(it));
        if(gen>=0) cout<<genId->at(gen)<<" - ";
        else cout<<"0 - ";
        cout<<trkPt->at(it)<<" "<<trkEta->at(it)<<" "<<trkPhi->at(it);
        if(gen>=0){cout<<"   M "; printMotherChain(gen);}
        cout<<endl;
    }

    cout<<endl<<endl<<"PV "<<iSsPV<<" "<<pvtNTracks->at(iSsPV)<<endl<<endl;
    for(auto it:trk_pv){
        if ( !(trkType->at( it ) & selFlag) ) continue; 
        if(!(( trkQuality->at( it ) >> 2 ) & 1)) continue;
        cout<<it<<" - ";
        int gen = GetClosestGen(trkEta->at(it), trkPhi->at(it), trkPt->at(it));
        if(gen>=0) cout<<genId->at(gen)<<" - ";
        else cout<<"0 - ";
        cout<<trkPt->at(it)<<" "<<trkEta->at(it)<<" "<<trkPhi->at(it);
        if(gen>=0){cout<<"   M "; printMotherChain(gen);}
        if(GetOverlappedTrack(it, &tkSsB) != -1) cout<<" <------- "<<endl;
        else cout<<endl;
        if(GetOverlappedTrack(it, &tkSsB) == -1)
            if(fabs(dXY(it, pvtX->at(iSsPV), pvtY->at(iSsPV)))>0.2){
                testList.push_back(it);
                hTest->Fill(dXY(it, pvtX->at(iSsPV), pvtY->at(iSsPV)));
            }
    }

    cout<<"Trk Charge "<<GetListCharge(&testList, 1.0)<<", lund "<<ssBLund<<endl;

    if(testList.size()>1){
        hmass_ssB_os->Fill(svtMass->at(iSsB), evtWeight);
        int isB = 0;
        if(GetListCharge(&testList, 1.0) < -QCut ) isB = +1;
        if(GetListCharge(&testList, 1.0) > +QCut ) isB = -1;

        if( TMath::Sign(1, ssBLund) == isB )
            hmass_ssB_osRT->Fill(svtMass->at(iSsB), evtWeight);

        if( TMath::Sign(1, ssBLund) == -1*isB )
            hmass_ssB_osWT->Fill(svtMass->at(iSsB), evtWeight);
    }
    return false;
*/
//-----------------------------------------JET-----------------------------------------

    int bestJet = -1;
    int bestJet_ = -1;
    float bestJetTag = CutDeepCSV;
    float bestJetPt = minPtJet;

    //SELECTION
    for (int iJet = 0; iJet<nJets; ++iJet){

        vector <int> jet_pfcs = pfCandFromJet( iJet );
        vector <int> jet_tks = tracksFromJet( iJet );

        //SELECTION
        if(goodJet(iJet)!=true) continue;
        if(jetNDau->at(iJet)<2) continue;
        if(abs(jetEta->at(iJet))>2.4) continue;
        if(jetPt->at(iJet)<minPtJet) continue;

        (tWriter->jetPt_v)->push_back(jetPt->at(iJet));
        (tWriter->jetCSV_v)->push_back(GetJetProbb(iJet));
        (tWriter->jetHasAncestor_v)->push_back(GetJetAncestor( iJet, &ListB ));

        float bTag = GetJetProbb(iJet);
        float cutbTag = CutDeepCSV;
        if(bTag < cutbTag) continue;

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
        if(abs(jetDz)>jetDzCut) continue;

        if(bTag>bestJetTag){
            bestJetTag = bTag;
            bestJet = iJet;
        }

        if(jetPt->at(iJet)>bestJetPt){
            bestJetPt = jetPt->at(iJet);
            bestJet_ = iJet;
        }

    }

    if(bestJetTag<0){
        bestJet = bestJet_;
    }

    //TAG
    if(bestJet < 0){
        (tWriter->osJet) = 0;
        (tWriter->evtNumber) = event_tot;
        if(saveNotTaggedEvets) tWriter->fill();
        return true;
    }

    (tWriter->osJet) = 1;

    hmass_ssB_os->Fill(svtMass->at(iSsB), evtWeight);

    //indices
    int iJet = bestJet;
    vector <int> jet_pfcs = pfCandFromJet( iJet );
    vector <int> jet_tks = tracksFromJet( iJet );

    //GENINFO
    int jetAncestor = GetJetAncestor( iJet, &ListB );

    //TAGGING VARIABLES
    //Separation
    cout<<endl;
    float jetDz=0;
    for(auto it:jet_tks) jetDz+= dZ(it, iSsPV);


    cout<<std::setprecision(3);
    cout<<" --- "<<iSsPV<<" ("<<pvtX->at(iSsPV)<<" "<<pvtY->at(iSsPV)<<" "<<pvtZ->at(iSsPV)<<")"<<endl;
    cout<<" --- "<<iJet<<" "<<jetPt->at(iJet)<<" "<<jetEta->at(iJet)<<" "<<jetPhi->at(iJet)<<"   "<<jetDz/jet_tks.size()<<endl;
    for(auto it:jet_tks){
        cout<<trkJet->at(it)<<" "<<trkPt->at(it)<<" "<<trkEta->at(it)<<" "<<trkPhi->at(it)<<" -  "<<dZ(it, iSsPV)<<" - ";
        cout<<trkPVtx->at(it)<<" ("<<pvtX->at(trkPVtx->at(it))<<" "<<pvtY->at(trkPVtx->at(it))<<" "<<pvtZ->at(trkPVtx->at(it))<<")"<<endl;
    }

    cout<<endl;

    float jetDzB = jetDz/jet_tks.size();
    float jetDrB = deltaR(jetEta->at( iJet ), jetPhi->at( iJet ), tB.Eta(), tB.Phi());

    //Jet Charge
    float jet_charge = GetJetCharge(iJet, kappa);

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

    if(jetAncestor<0)
        hmass_ssB_osJetwoAnc->Fill(svtMass->at(iSsB), evtWeight);
    else
        hmass_ssB_osJetwAnc->Fill(svtMass->at(iSsB), evtWeight);

    int isB = 0;
    (tWriter->osJetTag) = -2;

    if(jet_charge < -QCut ) isB = +1;
    if(jet_charge > +QCut ) isB = -1;

    if( TMath::Sign(1, ssBLund) == isB ){ 
        hmass_ssB_osRT->Fill(svtMass->at(iSsB), evtWeight);
        (tWriter->osJetTag) = 1;
    }

    if( TMath::Sign(1, ssBLund) == -1*isB ){
        hmass_ssB_osWT->Fill(svtMass->at(iSsB), evtWeight);
        (tWriter->osJetTag) = 0;
    }

    //------------------------------------------------TRACKS------------------------------------------------

    for(auto it:jet_tks){
        (tWriter->trkPt)->push_back(trkPt->at(it));
        (tWriter->trkEta)->push_back(trkEta->at(it));
        (tWriter->trkPhi)->push_back(trkPhi->at(it));
        (tWriter->trkCharge)->push_back(trkCharge->at(it));
        (tWriter->trkIsHighPurity)->push_back( ( trkQuality->at( it ) >> 2 ) & 1 );
        (tWriter->trkDxy)->push_back( dXYjet(it, iSsPV, iJet) );
        (tWriter->trkDz)->push_back(dZ(it, iSsPV));        
        (tWriter->trkIsInJet)->push_back( 1 );
    }

    for (int iTrk = 0; iTrk<nTracks; ++iTrk){

        if( deltaR(jetEta->at(iJet), jetPhi->at(iJet), trkEta->at(iTrk), trkPhi->at(iTrk)) > 0.5 ) continue;
        if(std::find(jet_tks.begin(), jet_tks.end(), iTrk) != jet_tks.end()) continue;

        (tWriter->trkPt)->push_back(trkPt->at(iTrk));
        (tWriter->trkEta)->push_back(trkEta->at(iTrk));
        (tWriter->trkPhi)->push_back(trkPhi->at(iTrk));
        (tWriter->trkCharge)->push_back(trkCharge->at(iTrk));
        (tWriter->trkIsHighPurity)->push_back( ( trkQuality->at( iTrk ) >> 2 ) & 1 );
        (tWriter->trkDxy)->push_back( dXYjet(iTrk, iSsPV, iJet) );
        (tWriter->trkDz)->push_back(dZ(iTrk, iSsPV));
        (tWriter->trkIsInJet)->push_back( 0 );        
    }

    (tWriter->evtNumber) = event_tot;
    tWriter->fill();

    return true;

}


void PDAnalyzer::endJob() {


// additional features

    tWriter->close();   // second ntuple

//printing results

    float eff   = static_cast<float>(hmass_ssB_os->Integral()) / static_cast<float>(hmass_ssB->Integral());
    float w     = static_cast<float>(hmass_ssB_osWT->Integral()) / static_cast<float>(hmass_ssB_os->Integral());
    float power = eff*pow(1-2*w, 2);

    cout<<"wAnc woAnc"<<endl;
    cout<<hmass_ssB_osJetwAnc->Integral()<<" "<<hmass_ssB_osJetwoAnc->Integral()<<endl<<endl;

    cout<<"#Bp eff% w% P%"<<endl;
    cout<<hmass_ssB->Integral()<<" "<<eff*100<<" "<<w*100<<" "<<power*100<<endl;

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
