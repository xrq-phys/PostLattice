#!/bin/bash

source VMC_Path.sh

for i in ./*.dir; do
	cd $i

	if [ ! -e output/zvo_out_101.dat ]; then

		#================================
		# Green Function Generation Phase
		cat > model.ini << EOF
[Control]
GreenOne = Green1.def
GreenTwo = Green2.def
Mode     = Gen
Verbose  = 0

[Physics]
System = Square
W      = 4

[Operator]
Spin_Structure = 1
AF_N_X         = 4
AF_N_Y         = 4
AF_Out         = sstruct.txt
SC             = d
SC_NumR        = 8
SC_Out         = sc.txt
EOF
		$anlys_ROOT/Source/Post model.ini

		#==================
		# Calculation Phase
		sed -i '.bak' '/.*InOrbital.*/d' namelist.def
		sed -i '.bak' 's/.*OneBodyG.*/OneBodyG Green1.def/' namelist.def
		sed -i '.bak' 's/.*TwoBodyG.*/TwoBodyG Green2.def/' namelist.def
		sed -i '.bak' 's/.*NDataIdxStart.*/NDataIdxStart 101/' modpara.def
		sed -i '.bak' 's/.*NVMCSample.*/NVMCSample 8000/' modpara.def
		sed -i '.bak' 's/.*NDataQtySmp.*/NDataQtySmp 8/' modpara.def
		sed -i '.bak' 's/.*NVMCCalMode.*/NVMCCalMode 1/' modpara.def
		$mVMC_ROOT/src/mVMC/vmc.out -e namelist.def output/zqp_opt.dat
	else
		echo "Already done GF generation for:"
		pwd
	fi

	if [ ! -e doublon.txt.101.dat ]; then

		#===============
		# Analysis Phase
		for j in output/zvo_cisajs_10*; do
			k=$(echo $j | sed -e "s/cisajs/cisajscktalt/") # Two body filename
			$anlys_ROOT/Source/Post model.ini $j $k
			k=$(echo $j | sed -e "s/output\/zvo_cisajs_/\./") # Number of file
			for l in ./*.txt; do
				mv $l "$l$k"
			done
		done
	else
		echo "Already done analysis for:"
		pwd
	fi

	cd ..
done

