export DSHOWER_ENV=/home/lthsmith/Desktop/dShowerOL
export DSHOWER_WITH_OL=1
export OPENLOOPS=/home/lthsmith/Desktop/dShowerOL/OpenLoops
export LD_LIBRARY_PATH=$OPENLOOPS/lib:$LD_LIBRARY_PATH

#python3 combinations.py

declare threads=14
declare proctrack=0
let "proctrac = 0"

make MC_MT_veto

file="/home/lthsmith/Desktop/dShowerOL/src/combs_WW.txt"
while IFS=',' read -r f1 f2 f3 f4 mult


do
	./MC_MT_veto 1 1000 $f1 $f2 $f3 $f4 0.0 &
	
	let "proctrac +=1"
	
	if [[ $((proctrac % threads)) -eq 0 ]]; then
		echo "The desired number of processes have been populated, holding until completed."
		wait
	fi
	
	
done<"$file"

#cd Integrals/Mixpol_plot
#python3 combiner_Mixpol.py 0.0
#cd ..
#cd ..
