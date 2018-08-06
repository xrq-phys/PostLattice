#!/bin/bash

source VMC_Path.sh

for i in 8 10 12 14 16; do
	mkdir U10N$i.dir
	cd    U10N$i.dir
	cat > model.def << EOF
W = 4
L = 4
model = "Fermion Hubbard"

lattice = "Square Lattice"
U = 10.0
t = 1.0
nelec = $i
EOF
	$mVMC_ROOT/src/mVMC/vmcdry.out model.def > init.log
	echo "Nelectron $(expr $i / 2)" >> modpara.def
	cd ..
done

