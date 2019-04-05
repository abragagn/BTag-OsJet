#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset() { autoReset(); }

PDSecondNtupleData() {

    trkPt           = new vector <float>;
    trkEta          = new vector <float>;
    trkPhi          = new vector <float>;
    trkCharge       = new vector <int>;
    trkDxy          = new vector <float>;
    trkDxySigned    = new vector <float>;
    trkDz           = new vector <float>;
    trkDxyz         = new vector <float>;
    trkExy          = new vector <float>;
    trkEz           = new vector <float>;
    trkExyz         = new vector <float>;
    trkDrJet        = new vector <float>;
    trkIsInJet      = new vector <int>;

}
virtual ~PDSecondNtupleData() {
}

void initTree() {
    treeName = "PDsecondTree";

    setBranch( "ssbMass", &ssbMass, "ssbMass/F", &b_ssbMass );
    setBranch( "ssbLund", &ssbLund, "ssbLund/I", &b_ssbLund );
    setBranch( "ssbIsTight", &ssbIsTight, "ssbIsTight/I", &b_ssbIsTight );

    setBranch( "ssbPVTx", &ssbPVTx, "ssbPVTx/F", &b_ssbPVTx );
    setBranch( "ssbPVTy", &ssbPVTy, "ssbPVTy/F", &b_ssbPVTy );
    setBranch( "ssbPVTz", &ssbPVTz, "ssbPVTz/F", &b_ssbPVTz );

    setBranch( "evtNumber", &evtNumber, "evtNumber/I", &b_evtNumber );
    setBranch( "evtWeight", &evtWeight, "evtWeight/I", &b_evtWeight );

    setBranch( "hltJpsiMu", &hltJpsiMu , "hltJpsiMu/I" , &b_hltJpsiMu );
    setBranch( "hltJpsiTrkTrk", &hltJpsiTrkTrk , "hltJpsiTrkTrk/I" , &b_hltJpsiTrkTrk );
    setBranch( "hltJpsiTrk", &hltJpsiTrk , "hltJpsiTrk/I" , &b_hltJpsiTrk );

    setBranch( "jetPt", &jetPt, "jetPt/F", &b_jetPt );
    setBranch( "jetEta", &jetEta, "jetEta/F", &b_jetEta );
    setBranch( "jetPhi", &jetPhi, "jetPhi/F", &b_jetPhi );
    setBranch( "jetCharge", &jetCharge, "jetCharge/F", &b_jetCharge );
    setBranch( "jetDeepCSV", &jetDeepCSV, "jetDeepCSV/F", &b_jetDeepCSV );
    setBranch( "jetDrB", &jetDrB, "jetDrB/F", &b_jetDrB );
    setBranch( "jetDzB", &jetDzB, "jetDzB/F", &b_jetDzB );

    // setBranch( "jetNDau", &jetNDau, "jetNDau/I", &b_jetNDau );
    // setBranch( "jetNHF", &jetNHF, "jetNHF/F", &b_jetNHF );
    // setBranch( "jetNEF", &jetNEF, "jetNEF/F", &b_jetNEF );
    // setBranch( "jetCHF", &jetCHF, "jetCHF/F", &b_jetCHF );
    // setBranch( "jetCEF", &jetCEF, "jetCEF/F", &b_jetCEF );
    // setBranch( "jetNCH", &jetNCH, "jetNCH/F", &b_jetNCH );

    setBranch( "jetHasAncestor", &jetHasAncestor, "jetHasAncestor/I", &b_jetHasAncestor );

    setBranch( "trkPt", &trkPt , 8192, 99, &b_trkPt );
    setBranch( "trkEta", &trkEta , 8192, 99, &b_trkEta );
    setBranch( "trkPhi", &trkPhi , 8192, 99, &b_trkPhi );
    setBranch( "trkCharge", &trkCharge , 8192, 99, &b_trkCharge );
    setBranch( "trkDxySigned", &trkDxySigned , 8192, 99, &b_trkDxySigned );
    setBranch( "trkDxy", &trkDxy , 8192, 99, &b_trkDxy );
    setBranch( "trkDz", &trkDz , 8192, 99, &b_trkDz );
    setBranch( "trkDxyz", &trkDxyz , 8192, 99, &b_trkDxyz );
    setBranch( "trkExy", &trkExy , 8192, 99, &b_trkExy );
    setBranch( "trkEz", &trkEz , 8192, 99, &b_trkEz );
    setBranch( "trkExyz", &trkExyz , 8192, 99, &b_trkExyz );
    setBranch( "trkDrJet", &trkDrJet , 8192, 99, &b_trkDrJet );
    setBranch( "trkIsInJet", &trkIsInJet , 8192, 99, &b_trkIsInJet );

}

float ssbMass, ssbPVTx, ssbPVTy, ssbPVTz;
int ssbLund, evtNumber, hltJpsiMu, hltJpsiTrkTrk, hltJpsiTrk, ssbIsTight, evtWeight;

TBranch *b_ssbMass, *b_ssbPVTx, *b_ssbPVTy, *b_ssbPVTz;
TBranch *b_ssbLund, *b_evtNumber, *b_hltJpsiMu, *b_hltJpsiTrkTrk, *b_hltJpsiTrk, *b_ssbIsTight, *b_evtWeight;

float jetPt, jetEta, jetPhi, jetCharge, jetCSV, jetDrB, jetDzB, jetDeepCSV;
float jetNHF, jetNEF, jetCHF, jetCEF, jetNCH;
int jetNDau, jetHasAncestor;

TBranch *b_jetPt, *b_jetEta, *b_jetPhi, *b_jetCharge, *b_jetCSV, *b_jetDrB, *b_jetDzB, *b_jetDeepCSV;
TBranch *b_jetNHF, *b_jetNEF, *b_jetCHF, *b_jetCEF, *b_jetNCH;
TBranch *b_jetNDau, *b_jetHasAncestor;

vector <float> *trkPt, *trkEta, *trkPhi;
vector <float> *trkDxy, *trkDxySigned, *trkDz, *trkDxyz, *trkExy, *trkEz, *trkExyz, *trkDrJet;
vector <int> *trkCharge, *trkIsInJet;
TBranch *b_trkPt, *b_trkEta, *b_trkPhi;
TBranch *b_trkDxy, *b_trkDxySigned, *b_trkDz, *b_trkDxyz, *b_trkExy, *b_trkEz, *b_trkExyz, *b_trkDrJet;
TBranch *b_trkCharge, *b_trkIsInJet;

private:

PDSecondNtupleData ( const PDSecondNtupleData& a );
PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

