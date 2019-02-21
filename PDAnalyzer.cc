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
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list hist.root -v outputFile ntu.root -v histoMode RECREATE -v use_gen t -n 10000
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

    getUserParameter( "saveNotTaggedEvets", saveNotTaggedEvets );
 
    getUserParameter( "minPtJet", minPtJet );
    getUserParameter( "kappa", kappa );
    getUserParameter( "CutCSV", CutCSV );
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

    int nc = 10;
    counter = new int[nc];
    for(int i=0;i<nc;i++)
        counter[i] = 0;

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
    hTest    = new TH1F( "hTest", "hTest", 100, 0, 30 );

    autoSavedObject =
    hTest2    = new TH1F( "hTest2", "hTest2", 100, 0, 30 );

    autoSavedObject =
    hTest3    = new TH1F( "hTest3", "hTest3", 100, -20, 20 );

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
    (tWriter->ssbPt) = tB.Pt();
    (tWriter->ssbEta) = tB.Eta();
    (tWriter->ssbPhi) = tB.Phi();
    (tWriter->ssbMass) = svtMass->at(ssbSVT);
    (tWriter->ssbIsTight) = isTight;

    (tWriter->ssbLxy) = GetCt2D(tB, ssbSVT) / (MassBs/tB.Pt());
    (tWriter->ssbCt2D) = GetCt2D(tB, ssbSVT);
    (tWriter->ssbCt2DErr) = GetCt2DErr(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt2DSigmaUnit) = GetCt2D(tB, ssbSVT, ssbPVT)/GetCt2DErr(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt3D) = GetCt3D(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt3DErr) = GetCt3DErr(tB, ssbSVT, ssbPVT);
    (tWriter->ssbCt3DSigmaUnit) = GetCt3D(tB, ssbSVT, ssbPVT)/GetCt3DErr(tB, ssbSVT, ssbPVT);

    (tWriter->ssbSVT) = ssbSVT;
    (tWriter->ssbPVT) = ssbPVT;
    
    (tWriter->ssbLund) = ssbLund;

    (tWriter->hltJpsiMu) = jpsimu;
    (tWriter->hltJpsiTrkTrk) = jpsitktk;
    (tWriter->hltJpsiTrk) = jpsitk;
    
    (tWriter->evtWeight) = evtWeight;
    (tWriter->evtNb) = ListLongLivedB.size();

    hmass_ssB->Fill(svtMass->at(ssbSVT), evtWeight);
    counter[0]+=evtWeight;

/*
    vector <int> trk_pv = tracksFromPV( ssbPVT );
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

    int bestos=-1;
    float bestospt=0;

    for( uint i=0; i<genId->size(); ++i ){
        if(!IsB(i)) continue;
        if(TagMixStatus(i)==2) continue;
        if((int)i==genBindex) continue;
        if(AreOverlapped( genPt->at(i), genEta->at(i), genPhi->at(i),  genPt->at(genBindex), genEta->at(genBindex), genPhi->at(genBindex) ) ) continue;
        if(deltaR(genEta->at(i), genPhi->at(i),genEta->at(genBindex), genPhi->at(genBindex))<=0.4) continue;
        if(genPt->at(i)>bestospt){
            bestospt = genPt->at(i);
            bestos = i;
        }
    }

    if(bestos==-1) return true;

    cout<<"------- PV "<<ssbPVT<<endl;

    const vector <int>& vD = allDaughters(bestos);
    cout<<genId->at(bestos)<<endl;
    for(uint i=0; i<vD.size(); ++i){
        cout<<genId->at(vD[i])<<" [ "<<genPt->at(vD[i])<<" "<<genEta->at(vD[i])<<" "<<genPhi->at(vD[i])<<" ]"<<endl;
    }
    cout<<endl;
    vector <int> trkCone;
    for(int i=0; i<nTracks; ++i){
        if(!(( trkQuality->at( i ) >> 2 ) & 1)) continue;
        if(deltaR(trkEta->at(i),trkPhi->at(i),genEta->at(bestos),genPhi->at(bestos))>0.4) continue;
        if(fabs(dZ(i, ssbPVT))>0.1) continue;
        trkCone.push_back(i);
        cout<<" [ "<<trkPt->at(i)<<" "<<trkEta->at(i)<<" "<<trkPhi->at(i)<<" ]  dz "<<dZ(i, ssbPVT)<<" ("<<trkPVtx->at(i)<<")"<<endl;
    }
    counter[0]++;

    if(trkCone.size()>=2){
        float q = GetListCharge(&trkCone, kappa);
        int d;
        if(q < -QCut ) d = +1;
        if(q > +QCut ) d = -1;

        if( TMath::Sign(1, ssbLund) == d )
            counter[1]++;

        if( TMath::Sign(1, ssbLund) == -1*d )
            counter[2]++;
    }

    return true;

/////////// DEBUG VERTICES

    vector<int> testList;

    cout<<endl<<"--------------------------------"<<endl;
    for( uint i=0; i<genId->size(); ++i ){
        if(!IsB(i)) continue;
        if(TagMixStatus(i)==2) continue;
        bool skip = false;
        for(auto it:testList)
            if(AreOverlapped( genPt->at(i), genEta->at(i), genPhi->at(i), 
                            genPt->at(it), genEta->at(it), genPhi->at(it) ) )
                skip = true;
        if(skip) continue;
        testList.push_back(i);
        cout<<std::setprecision(3);
        cout<<genId->at(i)<<" [ "<<genPt->at(i)<<" "<<genEta->at(i)<<" "<<genPhi->at(i)<<" ]"<<endl;
        printDaughterTreePt(i, " ");
        cout<<endl;
        int clJet = GetClosesJet(genPt->at(i), genEta->at(i), genPhi->at(i));
        if(clJet>=0){
            if(!AreOverlapped( genPt->at(i), genEta->at(i), genPhi->at(i),  genPt->at(genBindex), genEta->at(genBindex), genPhi->at(genBindex) ) ){
                hTest2->Fill(genPt->at(i));
            }
            cout<<"CLOSEST JET "<<clJet<<" [ "<<jetPt->at(clJet)<<" "<<jetEta->at(clJet)<<" "<<jetPhi->at(clJet)<<" ] dr "<<deltaR(jetEta->at(clJet), jetPhi->at(clJet), genEta->at(i), genPhi->at(i))<<endl;
            cout<<"CSV "<<GetJetProbb(clJet)<<endl;
            vector <int> jet_pfcs = pfCandFromJet( clJet );
            vector <int> jet_tks = tracksFromJet( clJet );
            for(auto it:jet_pfcs){
                int lund = GetClosestGen(pfcEta->at(it), pfcPhi->at(it), pfcPt->at(it));
                if(lund>=0) cout<<genId->at(lund);
                else cout<<"0";
                cout<<" ["<<pfcPt->at(it)<<" "<<pfcEta->at(it)<<" "<<pfcPhi->at(it)<<"] E "<<pfcE->at(it)<<", q "<<pfcCharge->at(it);
                if(pfcTrk->at(it)>=0) cout<<" dz "<<dZ(pfcTrk->at(it), ssbPVT)<<" ("<<trkPVtx->at(pfcTrk->at(it))<<")";
                cout<<endl;
            }
        }
        else if(!AreOverlapped( genPt->at(i), genEta->at(i), genPhi->at(i),  genPt->at(genBindex), genEta->at(genBindex), genPhi->at(genBindex) ) ) 
            hTest->Fill(genPt->at(i));

        cout<<endl<<endl;
    }

    int selFlag = PDEnumString::general | PDEnumString::packed;

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

    cout<<endl<<endl<<"PV "<<ssbPVT<<" "<<pvtNTracks->at(ssbPVT)<<endl<<endl;
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
            if(fabs(dXY(it, pvtX->at(ssbPVT), pvtY->at(ssbPVT)))>0.2){
                testList.push_back(it);
                hTest->Fill(dXY(it, pvtX->at(ssbPVT), pvtY->at(ssbPVT)));
            }
    }

    cout<<"Trk Charge "<<GetListCharge(&testList, 1.0)<<", lund "<<ssbLund<<endl;

    if(testList.size()>1){
        hmass_ssB_os->Fill(svtMass->at(ssbSVT), evtWeight);
        int isB = 0;
        if(GetListCharge(&testList, 1.0) < -QCut ) isB = +1;
        if(GetListCharge(&testList, 1.0) > +QCut ) isB = -1;

        if( TMath::Sign(1, ssbLund) == isB )
            hmass_ssB_osRT->Fill(svtMass->at(ssbSVT), evtWeight);

        if( TMath::Sign(1, ssbLund) == -1*isB )
            hmass_ssB_osWT->Fill(svtMass->at(ssbSVT), evtWeight);
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
        if(abs(jetEta->at(iJet))>2.5) continue;
        if(jetPt->at(iJet)<minPtJet) continue;

//        float bTag = GetJetProbb(iJet);
//        float cutbTag = CutDeepCSV;
        float bTag = jetCSV->at(iJet);
        float cutbTag = CutCSV;
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

    if(selectedJetTracks.size()==0){
        (tWriter->osJet) = 0;
        (tWriter->evtNumber) = event_tot;
        if(saveNotTaggedEvets) tWriter->fill();
        return true;
    }

    (tWriter->osJet) = 1;
    hmass_ssB_os->Fill(svtMass->at(ssbSVT), evtWeight);

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

//    cout<<"JET FOUND "<<" [ "<<jetPt->at(iJet)<<" "<<jetEta->at(iJet)<<" "<<jetPhi->at(iJet)<<" ]";
//    cout<<", dr "<<deltaR(jetEta->at(iJet), jetPhi->at(iJet), genEta->at(genBindex), genPhi->at(genBindex));
//    cout<<", q "<<jet_charge<<", lundSS "<<ssbLund<<endl;

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
        hmass_ssB_osJetwoAnc->Fill(svtMass->at(ssbSVT), evtWeight);
    else
        hmass_ssB_osJetwAnc->Fill(svtMass->at(ssbSVT), evtWeight);

    int isB = 0;
    (tWriter->osJetTag) = -2;

    if(jet_charge < -QCut ) isB = +1;
    if(jet_charge > +QCut ) isB = -1;

    if( TMath::Sign(1, ssbLund) == isB ){ 
        hmass_ssB_osRT->Fill(svtMass->at(ssbSVT), evtWeight);
        (tWriter->osJetTag) = 1;
        counter[1]+=evtWeight;
    }

    if( TMath::Sign(1, ssbLund) == -1*isB ){
        hmass_ssB_osWT->Fill(svtMass->at(ssbSVT), evtWeight);
        (tWriter->osJetTag) = 0;
        counter[2]+=evtWeight;
    }

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

    float eff   = static_cast<float>( counter[1]+counter[2] ) / static_cast<float>( counter[0] );
    float w     = static_cast<float>( counter[2] ) / static_cast<float>( counter[1]+counter[2] );
    float power = eff*pow(1-2*w, 2);

    cout<<counter[0]<<" - "<<counter[1]<<" - "<<counter[2]<<"     "<<std::setprecision(3)<<eff*100<<"     "<<w*100<<"     "<<power*100<<endl;

    return true;

}


void PDAnalyzer::endJob() {


// additional features

    tWriter->close();   // second ntuple

//printing results

    float eff   = static_cast<float>(hmass_ssB_os->Integral()) / static_cast<float>(hmass_ssB->Integral());
    float w     = static_cast<float>(hmass_ssB_osWT->Integral()) / static_cast<float>(hmass_ssB_os->Integral());
    float power = eff*pow(1-2*w, 2);

//    float eff   = static_cast<float>( counter[1]+counter[2] ) / static_cast<float>( counter[0] );
//    float w     = static_cast<float>( counter[2] ) / static_cast<float>( counter[1]+counter[2] );
//    float power = eff*pow(1-2*w, 2);

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
