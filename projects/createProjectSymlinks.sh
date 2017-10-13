
for FILE in "BPDE" "DamageSeq/0921" "DamageSeq/hiSeq" "DamageSeq/miseq" "DamageSeq/organs" "DamageSeq/ratio" "DamageSeq/UV" "ecoli/barcode" "ecoli/miseq" "ecoli/rnaSeq" "organsXR" "plant" "RNA-seq/organs" "XPA_ChipSeq" "XR-seq" "ZT"
do
    LINK="./$FILE/dataDir"
    rm -f $LINK
    ln -s /proj/sancarlb/users/ogun/$FILE $LINK
done