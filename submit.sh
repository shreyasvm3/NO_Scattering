#!/bin/bash

for i in {500..1000..500} 
do
	q=$(echo "scale = 3; $i/1000" | bc)
	mkdir -p Ei-$q/
	cd Ei-$q/
	if [ $i = 125 ]
	then 
		k=14
		t=$(echo "scale = 3; $k/1000" | bc)
		p=$((700000/$k))
	elif [ $i = 250 ] || [ $i = 375 ]
	then
		k=15
		t=$(echo "scale = 3; $k/1000" | bc)
		p=$((550000/$k))
	else
		k=16
		t=$(echo "scale = 3; $k/1000" | bc)
		p=$((400000/$k))

	fi
	for j in 12 13 14 15 17 18 
	do
		mkdir -p Nui-$j/
		cd Nui-$j/
		for n in 200
		do 
			mkdir -p Nb-$n/
			cd Nb-$n/
			for m in {4..4..1}
			do 
				o=$(($m+1))
				mkdir -p ntraj-10-$o/
				cd ntraj-10-$o/
		        	cp ../../../../dyn.x ../../../../run_multinode.sh ../../../../input .
		        	sed -i "s/input-ei/$q/g" input
		        	sed -i "s/input-ntraj/$m/g" input
		        	sed -i "s/input-nb/$n/g" input
		        	sed -i "s/input-nui/$j/g" input
		        	sed -i "s/input-dt/$t/g" input
		        	sed -i "s/input-nt/$p/g" input
				sbatch run_multinode.sh
				cd ../
			done
			cd ../
		done
		cd ../
	done
	cd ../
done
