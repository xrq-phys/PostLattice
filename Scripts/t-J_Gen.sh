#!/bin/bash

source VMC_Path.sh

for i in 14 16; do
	for j in 090 100 110 120 130 140; do
		for t in 0.5 1.0 1.5; do
			mkdir t$t-N$i-JRatio0.$j
			cd t$t-N$i-JRatio0.$j
			cat > model.inp << EOF
# NB: THIS IS NOT INPUT FILE FOR mVMC OR HPhi
# SEPARATE AND REMOVE COMMENT BEFORE RUNNING
&COMMON
 lattice = "Square Lattice"
 W       = 4
 L       = 4
 2Sz     = 0 # NOT NECESSARILY SINGLET!
&END
#
&HEIS
 model   = "Spin"
 J       = 1.0
 J'      = 0.1 # 0.1 J
&END
#
&HUB
 model   = "Fermion Hubbard"
 nelec   = 10
 t       = $t
&END
EOF
			# Concatinate Hamiltonian from Heisenberg and Hubbard model.
			$anlys_ROOT/t-j_input.py model.inp
			$mVMC_ROOT/src/mVMC/vmcdry.out model.heis.def > init.heis.log
			mv gutzwilleridx.def gutzwilleridx.def.bak
			$mVMC_ROOT/src/mVMC/vmcdry.out model.hubb.def > init.hubb.log
			mv gutzwilleridx.def.bak gutzwilleridx.def
			cat >> namelist.def << EOF
 Hund hund.def
 Exchange exchange.def
 InGutzwiller gutzwiller.def
EOF
			cat > gutzwiller.def << EOF
======================
NGutzwillerIdx   1
======================
==Excluding Doublons==
======================
    0   -10.0    0.0
EOF
			cd ..
		done
	done
done

