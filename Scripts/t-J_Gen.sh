#!/bin/bash

source VMC_Path.sh

$anlys_ROOT/t-j_input.py model.inp
$mVMC_ROOT/src/mVMC/vmcdry.out model.heis.def > init.heis.log
mv gutzwilleridx.def gutzwilleridx.def.bak
$mVMC_ROOT/src/mVMC/vmcdry.out model.hubb.def > init.hubb.log
mv gutzwilleridx.def.bak gutzwilleridx.def
cat >> namelist.def << EOF
 Hund  hund.def
 Exchange  exchange.def
EOF
