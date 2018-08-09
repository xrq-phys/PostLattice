#!/bin/bash

source VMC_Path.sh

for i in ./*.dir; do
	if [ ! -e $i/zqp_APOrbital_opt.dat ]; then
		cd $i
		echo "Working on $i"
		cat >> modpara.def << EOF
--------------------
Mix            0.4
EPS            10
IterationMax   20000
--------------------
EOF
		sed -i '.bak' '/.*NSRCG.*/d' modpara.def
		$mVMC_ROOT/src/ComplexUHF/UHF namelist.def
		if [ -e zqp_APOrbital_opt.dat ]; then
			echo ' InOrbital zqp_APOrbital_opt.dat' >> namelist.def
		fi
		sed -i '.bak' '/.*Mix.*/d' modpara.def
		sed -i '.bak' '/.*EPS.*/d' modpara.def
		sed -i '.bak' '/.*IterationMax.*/d' modpara.def
		echo ' NSRCG 0' >> modpara.def
		cd ..
	else
		echo "Already done for $i"
	fi
done

