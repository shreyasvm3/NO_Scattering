#!/bin/bash

for e in {125..1000..125}
do 
	f=$( echo "scale = 3; $e/1000" | bc)
	mkdir -p Ei-$f/
	cd Ei-$f/
        rm -r collected_data/
	mkdir -p collected_data/
	for l in 3 11 16
	do 
		for k in 150 175 200
		do
			cp Nui-$l/Nb-$k/ntraj-10-5/fort.1005 collected_data/Ei-$f-Nui-"$l"_RZ_Nb"$k".out 
			cp Nui-$l/Nb-$k/ntraj-10-5/fort.2005 collected_data/Ei-$f-Nui-"$l"_TP_Nb"$k".out 
			cp Nui-$l/Nb-$k/ntraj-10-5/fort.3005 collected_data/Ei-$f-Nui-"$l"_Pop_Nb"$k".out 
		done
	done
	cd collected_data/
	zip Ei-$f-collected-data.zip *
	cd ../
	cd ../
done
