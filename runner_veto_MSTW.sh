export DSHOWER_ENV=/home/lthsmith/Desktop/dShowerOL
export DSHOWER_WITH_OL=1
export OPENLOOPS=/home/lthsmith/Desktop/dShowerOL/OpenLoops
export LD_LIBRARY_PATH=$OPENLOOPS/lib:$LD_LIBRARY_PATH

#python3 combinations.py

declare threads=14
declare proctrack=0
let "proctrac = 0"
let "numfiles =0"

make MC_MT_veto

for i in $(seq 0.0 0.5 3.0);

do
	file="/home/lthsmith/Desktop/dShowerOL/src/combs_WW.txt"
	while IFS=',' read -r f1 f2 f3 f4 mult


	do	
		./MC_MT_veto 3 10000 $f1 $f2 $f3 $f4 $i &
		
		let "proctrac +=1"
		
		if [[ $((proctrac % threads)) -eq 0 ]]; then
			echo "The desired number of processes have been populated, holding until completed."
			wait
		fi
		
	
	done<"$file"
	
	let "numfiles +=1"
	
	cd Integrals/Temp_MSTW	
	python3 combiner_all.py $i
	cd ..
	cd ..
done

cd Integrals/Temp_MSTW
python3 Asym_plotter_MSTW.py $numfiles
cd ..
cd ..
