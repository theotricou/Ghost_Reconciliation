# cp ~/GitHub/Zombi/Parameters/* .




open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    printf '%.3d' $? >&3
    )&
}
N=6
open_sem $N


rm -rf sim_*
for transfer in 1; do
  for extinction in 0; do
    for nspecies in 20 40 60 80; do
      outdir=sim_$nspecies\_$extinction\_$transfer
      mkdir $outdir
      sed "s/TTTT/$transfer/g" GenomeParameters.tsv > $outdir/GenomeParameters.tsv
      sed "s/EEEE/$extinction/g; s/NNNN/$nspecies/g" SpeciesTreeParameters.tsv > $outdir/SpeciesTreeParameters.tsv
      cd $outdir
      for num in `seq 1 20`; do
        if [ ! -d  sim_$num ]; then
          bash /home/theo/GitHub/Ghost_Reconciliation/scripts/run_sim_rec.sh $num
        fi
      done
      cd ..
    done
  done
done


for transfer in 1; do
  for extinction in 0; do
    for nspecies in 100; do
      outdir=sim_$nspecies\_$extinction\_$transfer
      mkdir $outdir
      sed "s/TTTT/$transfer/g" GenomeParameters.tsv > $outdir/GenomeParameters.tsv
      sed "s/EEEE/$extinction/g; s/NNNN/$nspecies/g" SpeciesTreeParameters.tsv > $outdir/SpeciesTreeParameters.tsv
      cd $outdir
      for num in `seq 1 20`; do
        if [ ! -d  sim_$num ]; then
          bash /home/theo/GitHub/Ghost_Reconciliation/scripts/run_sim_rec.sh $num
        fi
      done
      cd ..
    done
  done
done





open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
     ( "$@"; )
    printf '%.3d' $? >&3
    )&
}
N=6
open_sem $N

for num in `seq 1 10`; do
  if [ ! -d  sim_$num ]; then
    run_with_lock bash /home/theo/GitHub/Ghost_Reconciliation/scripts/run_sim_rec.sh $num
  fi
done


seed=666
rep=100
rm -rf sim_*

python3 ~/GitHub/Zombi/Zombi.py T SpeciesTreeParameters.tsv sim_ $rep $seed 4
python3 ~/GitHub/Zombi/Zombi.py G GenomeParameters.tsv sim_ $rep $seed 4
python3 ~/GitHub/Zombi/SpeciesSampler.py n 20 sim_ $rep $seed 4



# GNU Ghost
