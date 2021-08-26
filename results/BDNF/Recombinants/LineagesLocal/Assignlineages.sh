#!/bin/bash

clear


# hyphy /home/aglucaci/hyphy-analyses/LabelTrees/label-tree.bf --tree /home/aglucaci/GeneInvestigator/results/BDNF/Recombinants/BDNF_codons_RDP_recombinationFree.fas.treefile --list /home/aglucaci/GeneInvestigator/results/BDNF/Recombinants/LineagesLocal/Artiodactyla.clade --output TEST.nwk --label TEST
WD=/home/aglucaci/GeneInvestigator/results/BDNF/Recombinants/LineagesLocal
TREE=/home/aglucaci/GeneInvestigator/results/BDNF/Recombinants/BDNF_codons_RDP_recombinationFree.fas.treefile
LabelTrees=/home/aglucaci/hyphy-analyses/LabelTrees/label-tree.bf


#CLADES
for clade in $WD/*.clade; do
   echo "#Clade species list: "$clade

   label="$(basename -- $clade)"
   label=${label%.clade}
   echo "#Label: "$label
 
   output="$TREE"."$label"
   echo "#Output: "$output

   echo ""

   echo hyphy $LabelTrees --tree $TREE --list $clade --output $output --label $label
   hyphy $LabelTrees --tree $TREE --list $clade --output $output --label $label
 
   echo ""
done
