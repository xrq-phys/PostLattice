#!/bin/bash

source VMC_Path.sh

for i in 16 14 12 10 8; do
	cd U10N$i
	cat >> modpara.def << EOF
--------------------
2Sz            0
Mix            0.5
EPS            20
IterationMax   1000
--------------------
EOF
	sed -i '.bak' '/.*NSRCG.*/d' modpara.def
	sed -i '.bak' '/.*Ncond.*/d' modpara.def
	$mVMC_ROOT/src/ComplexUHF/UHF namelist.def
	cat >> namelist.def << EOF
InOrbital zqp_APOrbital_opt.dat
EOF
	cat >> modpara.def << EOF
NSRCG 0
Ncond $i
EOF
	sed -i '.bak' '/.*Mix.*/d' modpara.def
	sed -i '.bak' '/.*EPS.*/d' modpara.def
	sed -i '.bak' '/.*IterationMax.*/d' modpara.def
	cd ..
done

