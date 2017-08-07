#~/usr/bin/bash

python mergeTSS.py
cd ../../XR-seq
python mergeTSS.py
cd -
cat dataDir/mergedTSS_NHF1_mergedRep.csv ../../XR-seq/dataDir/mergedTSS_NHF1_mergedRep.csv >dataDir/mergedTSS_NHF1_mergedRep_dsxr.csv
cat dataDir/mergedTES_NHF1_mergedRep.csv ../../XR-seq/dataDir/mergedTES_NHF1_mergedRep.csv >dataDir/mergedTES_NHF1_mergedRep_dsxr.csv