#!/bin/bash

#rm -f r1_reference_table7.csv  # Start fresh
#echo "transition_time,alpha,popsize,mat11,mat12,mat13,mat14,mat15,mat21,mat22,mat23,mat24,mat25,mat31,mat32,mat33,mat34,mat35,mat41,mat42,mat43,mat44,mat45,mat51,mat52,mat53,mat54,mat55,var_lastcolum,var_lastrow,eigen_var,zero_count" > r1_reference_table7.csv

for i in {1..76}; do
  echo "Simulation $i"

  python SIMgen.py

done &> r1_reference_table7_output.log &


echo "Script is running in the background. Check 'r1_reference_table7_output.log' for progress."

##nohup bash r1_referencetable7.sh &> r1_referencetable7_output.log &
##pkill -f r1_referencetable7.sh
