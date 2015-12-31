#!/bin/bash

rm -rf test;
rm -f results.txt

for mass in 500 900; do
  for fn in 0 1 2 9 12 19; do
    for var in   50  100  150  200  250  300  350  400  450  500  550  600  650  700  750  800  850  900  950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1450 1500 1550 1600 1650 1700 1750 1800 1850 1900 1950 2000; do
      echo 'TEST ' $mass $fn $var
      echo 'TEST ' $mass $fn $var >> results.txt
      cp wz_ana_test0.C wz_ana_test1.C
      export YYY=${var}
      sed -i s"/YYY/${YYY}/" wz_ana_test1.C
      source ~/EVAL_SH65 5_3_14
      rm -rf test/0;
      mkdir -p test/0;
      root -l -q -b wz_ana_test1.C+\($fn,4,\"ntuples_53x/backgroundEWK_skim14_sm_wz$mass.root\"\);
      mv histo_limits_qqwzll_shape_8TeV_Bin*.txt test/0;
      source ~/EVAL_SH66 6_1_1;
      combineCards.py -S test/0/histo_limits_qqwzll_shape_8TeV_Bin*.txt > test/0/qqwz.text;
      ~/releases/CMSSW_6_1_1/src/Combination/computeLimit.sh qqwz 0 exp $PWD/test | grep "Expected 5" >> results.txt
    done
  done
done
