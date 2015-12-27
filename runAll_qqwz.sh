#!/bin/bash

rm -rf cards_qqwz;

for mass in 500 900; do
  for fn in 0 1 2 9 12 19; do
    mkdir -p cards_qqwz/$mass/$fn;
    root -l -q -b wz_ana.C+\($fn,4,\"ntuples_53x/backgroundEWK_skim14_sm_wz$mass.root\"\);
    mv histo_limits_qqwzll_shape_8TeV_Bin*.txt cards_qqwz/$mass/$fn/;
  done
done
