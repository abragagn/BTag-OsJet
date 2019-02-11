#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h
#include <vector>
#include "NtuTool/Common/interface/TreeWrapper.h"
using namespace std;

class PDSecondNtupleData: public virtual TreeWrapper {

public:

void Reset() { autoReset(); }

PDSecondNtupleData() {


}
virtual ~PDSecondNtupleData() {
}

void initTree() {
    treeName = "PDsecondTree";

    setBranch( "ssbPt", &ssbPt, "ssbPt/F", &b_ssbPt );
    setBranch( "ssbEta", &ssbEta, "ssbEta/F", &b_ssbEta );
    setBranch( "ssbPhi", &ssbPhi, "ssbPhi/F", &b_ssbPhi );
    setBranch( "ssbDxy", &ssbDxy, "ssbDxy/F", &b_ssbDxy );
    setBranch( "ssbExy", &ssbExy, "ssbExy/F", &b_ssbExy );
    setBranch( "ssbDz", &ssbDz, "ssbDz/F", &b_ssbDz );
    setBranch( "ssbEz", &ssbEz, "ssbEz/F", &b_ssbEz );

    setBranch( "ssbMass", &ssbMass, "ssbMass/F", &b_ssbMass );
    setBranch( "jpsiMass", &jpsiMass, "jpsiMass/F", &b_jpsiMass );
    setBranch( "ssbSVT", &ssbSVT, "ssbSVT/I", &b_ssbSVT );
    setBranch( "ssbPVT", &ssbPVT, "ssbPVT/I", &b_ssbPVT );
    setBranch( "ssbLund", &ssbLund, "ssbLund/I", &b_ssbLund );
    setBranch( "ssbIsTight", &ssbIsTight, "ssbIsTight/I", &b_ssbIsTight );

    setBranch( "ssbLxy", &ssbLxy , "ssbLxy/F" , &b_ssbLxy );
    setBranch( "ssbCt2D", &ssbCt2D , "ssbCt2D/F" , &b_ssbCt2D );
    setBranch( "ssbCt2DErr", &ssbCt2DErr , "ssbCt2DErr/F" , &b_ssbCt2DErr );
    setBranch( "ssbCt2DSigmaUnit", &ssbCt2DSigmaUnit , "ssbCt2DSigmaUnit/F" , &b_ssbCt2DSigmaUnit );
    setBranch( "ssbCt3D", &ssbCt3D , "ssbCt3D/F" , &b_ssbCt3D );
    setBranch( "ssbCt3DErr", &ssbCt3DErr , "ssbCt3DErr/F" , &b_ssbCt3DErr );
    setBranch( "ssbCt3DSigmaUnit", &ssbCt3DSigmaUnit , "ssbCt3DSigmaUnit/F" , &b_ssbCt3DSigmaUnit );

    setBranch( "evtNumber", &evtNumber, "evtNumber/I", &b_evtNumber );
    setBranch( "evtWeight", &evtWeight, "evtWeight/I", &b_evtWeight );
    setBranch( "evtNb", &evtNb, "evtNb/I", &b_evtNb );

    setBranch( "hltJpsiMu", &hltJpsiMu , "hltJpsiMu/I" , &b_hltJpsiMu );
    setBranch( "hltJpsiTrkTrk", &hltJpsiTrkTrk , "hltJpsiTrkTrk/I" , &b_hltJpsiTrkTrk );
    setBranch( "hltJpsiTrk", &hltJpsiTrk , "hltJpsiTrk/I" , &b_hltJpsiTrk );

    setBranch( "osJetTag", &osJetTag, "osJetTag/I", &b_osJetTag );      
    setBranch( "osJet", &osJet, "osJet/I", &b_osJet );      

    setBranch( "jetPt", &jetPt, "jetPt/F", &b_jetPt );
    setBranch( "jetEta", &jetEta, "jetEta/F", &b_jetEta );
    setBranch( "jetPhi", &jetPhi, "jetPhi/F", &b_jetPhi );
    setBranch( "jetCharge", &jetCharge, "jetCharge/F", &b_jetCharge );
    setBranch( "jetCSV", &jetCSV, "jetCSV/F", &b_jetCSV );
    setBranch( "jetDFprobb", &jetDFprobb, "jetDFprobb/F", &b_jetDFprobb );
    setBranch( "jetDrB", &jetDrB, "jetDrB/F", &b_jetDrB );
    setBranch( "jetDzB", &jetDzB, "jetDzB/F", &b_jetDzB );

    setBranch( "jetSize", &jetSize, "jetSize/I", &b_jetSize );
    setBranch( "jetHasAncestor", &jetHasAncestor, "jetHasAncestor/I", &b_jetHasAncestor );

}

float ssbPt, ssbEta, ssbPhi, ssbMass, jpsiMass, ssbDxy, ssbExy, ssbDz, ssbEz;
float ssbLxy, ssbCt2D, ssbCt2DErr, ssbCt2DSigmaUnit, ssbCt3D, ssbCt3DErr, ssbCt3DSigmaUnit;
int ssbSVT, ssbPVT, ssbLund, evtNumber, hltJpsiMu, hltJpsiTrkTrk, hltJpsiTrk, ssbIsTight, evtNb, evtWeight;

TBranch *b_ssbPt, *b_ssbEta, *b_ssbPhi, *b_ssbMass, *b_jpsiMass, *b_ssbDxy, *b_ssbExy, *b_ssbDz, *b_ssbEz;
TBranch *b_ssbLxy, *b_ssbCt2D, *b_ssbCt2DErr, *b_ssbCt2DSigmaUnit, *b_ssbCt3D, *b_ssbCt3DErr, *b_ssbCt3DSigmaUnit;
TBranch *b_ssbSVT, *b_ssbPVT, *b_ssbLund, *b_evtNumber, *b_hltJpsiMu, *b_hltJpsiTrkTrk, *b_hltJpsiTrk, *b_ssbIsTight, *b_evtWeight, *b_evtNb;

int osJetTag, osJet;
TBranch *b_osJetTag, *b_osJet;

float jetPt, jetEta, jetPhi, jetCharge, jetCSV, jetDrB, jetDzB, jetDFprobb;
int jetSize, jetHasAncestor;

TBranch *b_jetPt, *b_jetEta, *b_jetPhi, *b_jetCharge, *b_jetCSV, *b_jetDrB, *b_jetDzB, *b_jetDFprobb;
TBranch *b_jetSize, *b_jetHasAncestor;

private:

PDSecondNtupleData ( const PDSecondNtupleData& a );
PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

