#!/bin/bash

source VMC_Path.sh

for i in ./U10*; do
	cd $i

	if [ ! -e output/zvo_out_201.dat ]; then

		#================================
		# Green Function Generation Phase
		cat > GreenSC.def << EOF
=======================
 NCisAjsCktAlt   XXXXX
=======================
======== SC ===========
=======================
EOF
		$anlys_ROOT/opt_generators/sc_green.cc.x 4 4 >> GreenSC.def
		nGRNPT=$(awk '/./{line=$0} END{print line}' GreenSC.def)
		sed -i '.bak' "/.*NCisAjsCktAltSC.*/d" GreenSC.def
		sed -i '.bak' "s/.*NCisAjsCktAlt.*/$nGRNPT/" GreenSC.def

		#==================
		# Calculation Phase
		sed -i '.bak' 's/.*NVMCCalMode.*/NVMCCalMode 1/' modpara.def
		sed -i '.bak' 's/.*NDataQtySmp.*/NDataQtySmp 8/' modpara.def
		sed -i '.bak' 's/.*NVMCSample.*/NVMCSample 8000/' modpara.def
		sed -i '.bak' 's/.*NDataIdxStart.*/NDataIdxStart 201/' modpara.def
		sed -i '.bak' 's/.*TwoBodyG.*/TwoBodyG GreenSC.def/' namelist.def
		sed -i '.bak' '/.*InOrbital.*/d' namelist.def
		$mVMC_ROOT/src/mVMC/vmc.out -e namelist.def output/zqp_opt.dat

	else
		echo "Already done GF generation for:"
		pwd
	fi

	if [ ! -e sc.txt.001.dat.sc ]; then

		#===============
		# Analysis Phase
		latt_w=4
		latt_r=8
		for j in output/zvo_cisajs_20*; do
			k=$(echo $j | sed -e "s/cisajs/cisajscktalt/") # Two body filename
			$anlys_ROOT/opt_generators/sc_measure.cc.x $latt_w $j $k $latt_r
			k=$(echo $j | sed -e "s/output\/zvo_cisajs_/\./") # Number of file
			for l in ./*.txt; do
				mv $l "$l$k.sc"
			done
		done
	else
		echo "Already done SC analysis for:"
		pwd
	fi

	cd ..
done

