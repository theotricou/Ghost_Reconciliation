#!/bin/bash

sim=$1

python3 ~/GitHub/Zombi/Zombi.py T SpeciesTreeParameters.tsv sim_$sim

python3 ~/GitHub/Zombi/Zombi.py G GenomeParameters.tsv sim_$sim

cd sim_$sim

grep -P "\tT" G/Gene_families/* | awk '{print $1, $3}' |
  sed "s/.tsv.* / /g" | sed "s/;/ /g" | awk '{print $1, $2, $6}' |
  sed "s@.*/@@g" | sed 's/events/prunedtree/g' > Transfers.txt

NUM=1
while [ ! -f aletree ]; do
  if grep -q "(" G/Gene_trees/"${NUM}"_prunedtree.nwk;
  then
    ALEobserve G/Gene_trees/"${NUM}"_prunedtree.nwk
    ALEml_undated T/ExtantTree.nwk G/Gene_trees/"${NUM}"_prunedtree.nwk.ale output_species_tree=y sample=0
    mv ExtantTree.nwk_"${NUM}"_prunedtree.nwk.ale.spTree aletree
    rm ExtantTree.nwk*
  else
   let NUM++
  fi
done

for i in G/Gene_trees/*_prunedtree.nwk; do
  if grep -q "(" "$i"; then
    gene=`echo "${i/*\/}"`
    # ALE reconciliation
    if [ ! -f $gene.ale.uTs ]; then
      ALEobserve $i
      ALEml_undated T/ExtantTree.nwk $i.ale
      sed 's/([^)]*)//g' ExtantTree.nwk_$gene.ale.uTs | grep -o "[0-9].*" > $gene.ale.uTs
      rm ExtantTree.nwk_$gene.ale.*
    fi
    # RANGER reconciliation
    if [ ! -f "${gene}".ranger.uTs ]; then
      # Create a file with the species tree and the gene tree
      cat aletree "${i/.ale}" > "${gene}"_ranger_input
      ~/Downloads/9_Tools/RANGER-DTL/CorePrograms/Ranger-DTL.linux -i "${gene}"_ranger_input -o "${gene}"_ranger_output
      grep "Transfer," "${gene}"_ranger_output | awk '{print $8, $11}' | sed "s/,//g" | sed "s/$/ 1/g" > "${gene}".ranger.uTs
      rm "${gene}"_ranger*
    fi
    # Notung reconciliation
    if [ ! -f ${gene}.notung.uTs ]; then
      # Transform distances from gene tree in bootstraps for rooting error when using --parsable
      # sed "s/:[0-9]\.[0-9]*/:100/g; s/:[0-9]*\;/:100\;/g" $i > ${gene}_notung_input
      java -jar ~/Downloads/9_Tools/Notung-2.9.1.5/Notung-2.9.1.5.jar --reconcile -g $i -s aletree --bootstraps name --infertransfers true --multsols 100 --parsable
      grep "#T" ${gene}.reconciled.0.parsable.txt | awk '{print $4,$5}' | grep -o "[0-9].*" | sed "s/$/ 1/g" > ${gene}.notung.uTs
      rm ${gene}.reconciled.0*
    fi
  fi
done

# remove empty files
find . -type f -empty -print -delete

# cat all in one
for i in *uTs; do
  gene=`cut -d'.' -f 1 <<< $i`
  rec=`cut -d'.' -f 3 <<< $i`
  sed "s/^/$gene $rec /g" $i
done > Transfers_rec.txt

python3 /home/theo/GitHub/Ghost_Reconciliation/scripts/nodes_mappings.py

# GNU Ghost
