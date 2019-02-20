#!/bin/sh
#BSUB -o testHevjin.log
eval `scram runtime -sh`
pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BsToJpsiPhi_2017_DCAP.list histHevjin.root -v outputFile ntuHevjin.root -v histoMode RECREATE -v use_gen t -v useHLT f -n 10000000
