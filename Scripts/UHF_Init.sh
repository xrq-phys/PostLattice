#!/bin/bash

source VMC_Path.sh

for i in ./*.dir; do
	cd $i
	cat >> modpara.def << EOF
--------------------
Mix            0.4
EPS            10
IterationMax   8000
--------------------
EOF
	sed -i '.bak' '/.*NSRCG.*/d' modpara.def
	$mVMC_ROOT/src/ComplexUHF/UHF namelist.def
	echo ' InOrbital zqp_APOrbital_opt.dat' >> namelist.def
	sed -i '.bak' '/.*Mix.*/d' modpara.def
	sed -i '.bak' '/.*EPS.*/d' modpara.def
	sed -i '.bak' '/.*IterationMax.*/d' modpara.def
	echo ' NSRCG 0' >> modpara.def
	cd ..
done

