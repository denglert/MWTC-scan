#!/bin/bash

TAG=$1
CARD=./Cards/param_card_${TAG}.dat

# /* TCInput Parameters */
gt=2
MA=500
PS=0.3
rs=0.
MH=200.


./MWT_Calculator ${gt} ${MA} ${PS} ${rs} ${MH} > ${CARD}
cat ${CARD}
