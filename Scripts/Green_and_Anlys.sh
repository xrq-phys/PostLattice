#!/bin/bash

source VMC_Path.sh

for i in ./U10*; do
	cd $i

	if [ ! -e output/zvo_out_101.dat ]; then

		#================================
		# Green Function Generation Phase
		cat > GreenSC.def << EOF
=======================
 NCisAjsCktAlt   XXXXX
=======================
======== SC ===========
=======================
EOF

		$anlys_ROOT/Tools_CC/sc_green.cc.x 4 4 >> GreenSC.def
		nGRNPT=$(awk '/./{line=$0} END{print line}' GreenSC.def)
		sed -i '.bak' "/.*NCisAjsCktAltSC.*/d" GreenSC.def
		sed -i '.bak' "s/.*NCisAjsCktAlt.*/$nGRNPT/" GreenSC.def
		$anlys_ROOT/Scripts/Merge_Green.sh greentwo.def GreenSC.def GreenFull.def

		#==================
		# Calculation Phase
		sed -i '.bak' '/.*InOrbital.*/d' namelist.def
		sed -i '.bak' 's/.*TwoBodyG.*/TwoBodyG GreenFull.def/' namelist.def
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
		cat > system.py << EOF
sys_type = "square2d"
a_length = [4, 4]
runmode = "wick"
flip_index = []
for i in range(16):
    flip_index.extend([ 0 ] * 32)
    flip_index.extend([ 1 ] * 32)
    flip_index.extend([ 0 ] * 32)
EOF
		latt_w=4
		latt_r=8
		for j in output/zvo_cisajs_10*; do
			k=$(echo $j | sed -e "s/cisajs/cisajscktalt/") # Two body filename
			$anlys_ROOT/tokyo_face.py system.py $j $k -
			$anlys_ROOT/Tools_CC/sc_measure.cc.x $latt_w $j $k $latt_r
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

