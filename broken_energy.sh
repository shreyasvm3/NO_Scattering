#!/bin/bash

echo "Nui = 3", "Nui = 11", "Nui = 16" >> energy_stats.txt
echo "Nui = 3", "Nui = 11", "Nui = 16" >> invpot_stats.txt
for e in {125..1000..125}
do 
	f=$( echo "scale = 3; $e/1000" | bc)
	aa=$(grep "energy" Ei-$f/Nui-3/Nb-200/ntraj-10-5/Traj_Info.out | tail -1 | awk '{split($0,a,":"); print a[2]}') 
        bb=$(grep "energy" Ei-$f/Nui-3/Nb-200/ntraj-10-5/Traj_Info.out | tail -1 | awk '{split($0,a,":"); print a[2]}') 
        cc=$(grep "energy" Ei-$f/Nui-3/Nb-200/ntraj-10-5/Traj_Info.out | tail -1 | awk '{split($0,a,":"); print a[2]}') 
	aaa=$(grep "inv" Ei-$f/Nui-3/Nb-200/ntraj-10-5/Traj_Info.out | tail -1 | awk '{split($0,a,":"); print a[2]}') 
        bbb=$(grep "inv" Ei-$f/Nui-3/Nb-200/ntraj-10-5/Traj_Info.out | tail -1 | awk '{split($0,a,":"); print a[2]}') 
        ccc=$(grep "inv" Ei-$f/Nui-3/Nb-200/ntraj-10-5/Traj_Info.out | tail -1 | awk '{split($0,a,":"); print a[2]}') 
	echo "Ei-$f", "$aa", "$bb", "$cc" >> energy_stats.txt
	echo "Ei-$f", "$aaa", "$bbb", "$ccc" >> invpot_stats.txt
done
