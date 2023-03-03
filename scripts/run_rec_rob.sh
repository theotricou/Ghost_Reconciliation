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

for num in `seq 1 10`; do
  if [ ! -d  sim_$num ]; then
    run_with_lock bash /home/theo/GitHub/Ghost_Reconciliation/scripts/run_sim_rec.sh $num
  fi
done
