#!/bin/bash

#rm -f r1_reference_table75.csv  # Start fresh
#echo "sim_number,transition_time,alpha,popsize,mean_psi,var_psi,mean_psi2,var_psi2,mean_psi3,var_psi3,mean_B1,var_B1,mean_B2,var_B2,mean_Tpi,var_Tpi,mean_Wtheta,var_Wtheta,mean_TajimasD,var_TajimasD,Delta1,Delta2,Delta3,Delta4,Delta5,Delta6,Delta7,Delta8,Delta9,r_EMIBD,var_Delta1,var_Delta2,var_Delta3,var_Delta4,var_Delta5,var_Delta6,var_Delta7,var_Delta8,var_Delta9,var_r_EMIBD,F_Mean,F_SD,var_F_Mean,var_F_SD,var2_Delta1,var2_Delta2,var2_Delta3,var2_Delta4,var2_Delta5,var2_Delta6,var2_Delta7,var2_Delta8,var2_Delta9,var2_r_EMIBD,var2_F_Mean,var2_F_SD,mat11_Wtheta,mat11_psi,mat11_psi2,mat11_Delta1,mat11_Delta2,mat11_Delta3,mat11_Delta4,mat11_Delta5,mat11_Delta6,mat11_Delta7,mat11_Delta8,mat11_Delta9,mat11_r_EMIBD,mat11_F_Mean,mat11_F_SD,mat12_Wtheta,mat12_psi,mat12_psi2,mat12_Delta1,mat12_Delta2,mat12_Delta3,mat12_Delta4,mat12_Delta5,mat12_Delta6,mat12_Delta7,mat12_Delta8,mat12_Delta9,mat12_r_EMIBD,mat12_F_Mean,mat12_F_SD,mat13_Wtheta,mat13_psi,mat13_psi2,mat13_Delta1,mat13_Delta2,mat13_Delta3,mat13_Delta4,mat13_Delta5,mat13_Delta6,mat13_Delta7,mat13_Delta8,mat13_Delta9,mat13_r_EMIBD,mat13_F_Mean,mat13_F_SD,mat14_Wtheta,mat14_psi,mat14_psi2,mat14_Delta1,mat14_Delta2,mat14_Delta3,mat14_Delta4,mat14_Delta5,mat14_Delta6,mat14_Delta7,mat14_Delta8,mat14_Delta9,mat14_r_EMIBD,mat14_F_Mean,mat14_F_SD,mat15_Wtheta,mat15_psi,mat15_psi2,mat15_Delta1,mat15_Delta2,mat15_Delta3,mat15_Delta4,mat15_Delta5,mat15_Delta6,mat15_Delta7,mat15_Delta8,mat15_Delta9,mat15_r_EMIBD,mat15_F_Mean,mat15_F_SD,mat21_Wtheta,mat21_psi,mat21_psi2,mat21_Delta1,mat21_Delta2,mat21_Delta3,mat21_Delta4,mat21_Delta5,mat21_Delta6,mat21_Delta7,mat21_Delta8,mat21_Delta9,mat21_r_EMIBD,mat21_F_Mean,mat21_F_SD,mat22_Wtheta,mat22_psi,mat22_psi2,mat22_Delta1,mat22_Delta2,mat22_Delta3,mat22_Delta4,mat22_Delta5,mat22_Delta6,mat22_Delta7,mat22_Delta8,mat22_Delta9,mat22_r_EMIBD,mat22_F_Mean,mat22_F_SD,mat23_Wtheta,mat23_psi,mat23_psi2,mat23_Delta1,mat23_Delta2,mat23_Delta3,mat23_Delta4,mat23_Delta5,mat23_Delta6,mat23_Delta7,mat23_Delta8,mat23_Delta9,mat23_r_EMIBD,mat23_F_Mean,mat23_F_SD,mat24_Wtheta,mat24_psi,mat24_psi2,mat24_Delta1,mat24_Delta2,mat24_Delta3,mat24_Delta4,mat24_Delta5,mat24_Delta6,mat24_Delta7,mat24_Delta8,mat24_Delta9,mat24_r_EMIBD,mat24_F_Mean,mat24_F_SD,mat25_Wtheta,mat25_psi,mat25_psi2,mat25_Delta1,mat25_Delta2,mat25_Delta3,mat25_Delta4,mat25_Delta5,mat25_Delta6,mat25_Delta7,mat25_Delta8,mat25_Delta9,mat25_r_EMIBD,mat25_F_Mean,mat25_F_SD,mat31_Wtheta,mat31_psi,mat31_psi2,mat31_Delta1,mat31_Delta2,mat31_Delta3,mat31_Delta4,mat31_Delta5,mat31_Delta6,mat31_Delta7,mat31_Delta8,mat31_Delta9,mat31_r_EMIBD,mat31_F_Mean,mat31_F_SD,mat32_Wtheta,mat32_psi,mat32_psi2,mat32_Delta1,mat32_Delta2,mat32_Delta3,mat32_Delta4,mat32_Delta5,mat32_Delta6,mat32_Delta7,mat32_Delta8,mat32_Delta9,mat32_r_EMIBD,mat32_F_Mean,mat32_F_SD,mat33_Wtheta,mat33_psi,mat33_psi2,mat33_Delta1,mat33_Delta2,mat33_Delta3,mat33_Delta4,mat33_Delta5,mat33_Delta6,mat33_Delta7,mat33_Delta8,mat33_Delta9,mat33_r_EMIBD,mat33_F_Mean,mat33_F_SD,mat34_Wtheta,mat34_psi,mat34_psi2,mat34_Delta1,mat34_Delta2,mat34_Delta3,mat34_Delta4,mat34_Delta5,mat34_Delta6,mat34_Delta7,mat34_Delta8,mat34_Delta9,mat34_r_EMIBD,mat34_F_Mean,mat34_F_SD,mat35_Wtheta,mat35_psi,mat35_psi2,mat35_Delta1,mat35_Delta2,mat35_Delta3,mat35_Delta4,mat35_Delta5,mat35_Delta6,mat35_Delta7,mat35_Delta8,mat35_Delta9,mat35_r_EMIBD,mat35_F_Mean,mat35_F_SD,mat41_Wtheta,mat41_psi,mat41_psi2,mat41_Delta1,mat41_Delta2,mat41_Delta3,mat41_Delta4,mat41_Delta5,mat41_Delta6,mat41_Delta7,mat41_Delta8,mat41_Delta9,mat41_r_EMIBD,mat41_F_Mean,mat41_F_SD,mat42_Wtheta,mat42_psi,mat42_psi2,mat42_Delta1,mat42_Delta2,mat42_Delta3,mat42_Delta4,mat42_Delta5,mat42_Delta6,mat42_Delta7,mat42_Delta8,mat42_Delta9,mat42_r_EMIBD,mat42_F_Mean,mat42_F_SD,mat43_Wtheta,mat43_psi,mat43_psi2,mat43_Delta1,mat43_Delta2,mat43_Delta3,mat43_Delta4,mat43_Delta5,mat43_Delta6,mat43_Delta7,mat43_Delta8,mat43_Delta9,mat43_r_EMIBD,mat43_F_Mean,mat43_F_SD,mat44_Wtheta,mat44_psi,mat44_psi2,mat44_Delta1,mat44_Delta2,mat44_Delta3,mat44_Delta4,mat44_Delta5,mat44_Delta6,mat44_Delta7,mat44_Delta8,mat44_Delta9,mat44_r_EMIBD,mat44_F_Mean,mat44_F_SD,mat45_Wtheta,mat45_psi,mat45_psi2,mat45_Delta1,mat45_Delta2,mat45_Delta3,mat45_Delta4,mat45_Delta5,mat45_Delta6,mat45_Delta7,mat45_Delta8,mat45_Delta9,mat45_r_EMIBD,mat45_F_Mean,mat45_F_SD,mat51_Wtheta,mat51_psi,mat51_psi2,mat51_Delta1,mat51_Delta2,mat51_Delta3,mat51_Delta4,mat51_Delta5,mat51_Delta6,mat51_Delta7,mat51_Delta8,mat51_Delta9,mat51_r_EMIBD,mat51_F_Mean,mat51_F_SD,mat52_Wtheta,mat52_psi,mat52_psi2,mat52_Delta1,mat52_Delta2,mat52_Delta3,mat52_Delta4,mat52_Delta5,mat52_Delta6,mat52_Delta7,mat52_Delta8,mat52_Delta9,mat52_r_EMIBD,mat52_F_Mean,mat52_F_SD,mat53_Wtheta,mat53_psi,mat53_psi2,mat53_Delta1,mat53_Delta2,mat53_Delta3,mat53_Delta4,mat53_Delta5,mat53_Delta6,mat53_Delta7,mat53_Delta8,mat53_Delta9,mat53_r_EMIBD,mat53_F_Mean,mat53_F_SD,mat54_Wtheta,mat54_psi,mat54_psi2,mat54_Delta1,mat54_Delta2,mat54_Delta3,mat54_Delta4,mat54_Delta5,mat54_Delta6,mat54_Delta7,mat54_Delta8,mat54_Delta9,mat54_r_EMIBD,mat54_F_Mean,mat54_F_SD,mat55_Wtheta,mat55_psi,mat55_psi2,mat55_Delta1,mat55_Delta2,mat55_Delta3,mat55_Delta4,mat55_Delta5,mat55_Delta6,mat55_Delta7,mat55_Delta8,mat55_Delta9,mat55_r_EMIBD,mat55_F_Mean,mat55_F_SD,var_lastcolumn_Wtheta,var_lastcolumn_psi,var_lastcolumn_psi2,var_lastcolumn_Delta1,var_lastcolumn_Delta2,var_lastcolumn_Delta3,var_lastcolumn_Delta4,var_lastcolumn_Delta5,var_lastcolumn_Delta6,var_lastcolumn_Delta7,var_lastcolumn_Delta8,var_lastcolumn_Delta9,var_lastcolumn_r_EMIBD,var_lastcolumn_F_Mean,var_lastcolumn_F_SD,var_lastrow_Wtheta,var_lastrow_psi,var_lastrow_psi2,var_lastrow_Delta1,var_lastrow_Delta2,var_lastrow_Delta3,var_lastrow_Delta4,var_lastrow_Delta5,var_lastrow_Delta6,var_lastrow_Delta7,var_lastrow_Delta8,var_lastrow_Delta9,var_lastrow_r_EMIBD,var_lastrow_F_Mean,var_lastrow_F_SD,eigen_var_Wtheta,eigen_var_psi,eigen_var_psi2,eigen_var_Delta1,eigen_var_Delta2,eigen_var_Delta3,eigen_var_Delta4,eigen_var_Delta5,eigen_var_Delta6,eigen_var_Delta7,eigen_var_Delta8,eigen_var_Delta9,eigen_var_r_EMIBD,eigen_var_F_Mean,eigen_var_F_SD,zero_count_Wtheta,zero_count_psi,zero_count_psi2,zero_count_Delta1,zero_count_Delta2,zero_count_Delta3,zero_count_Delta4,zero_count_Delta5,zero_count_Delta6,zero_count_Delta7,zero_count_Delta8,zero_count_Delta9,zero_count_r_EMIBD,zero_count_F_Mean,zero_count_F_SD,ht1_Wtheta,ht1_psi,ht1_psi2,ht1_Delta1,ht1_Delta2,ht1_Delta3,ht1_Delta4,ht1_Delta5,ht1_Delta6,ht1_Delta7,ht1_Delta8,ht1_Delta9,ht1_r_EMIBD,ht1_F_Mean,ht1_F_SD,ht2_Wtheta,ht2_psi,ht2_psi2,ht2_Delta1,ht2_Delta2,ht2_Delta3,ht2_Delta4,ht2_Delta5,ht2_Delta6,ht2_Delta7,ht2_Delta8,ht2_Delta9,ht2_r_EMIBD,ht2_F_Mean,ht2_F_SD,ht3_Wtheta,ht3_psi,ht3_psi2,ht3_Delta1,ht3_Delta2,ht3_Delta3,ht3_Delta4,ht3_Delta5,ht3_Delta6,ht3_Delta7,ht3_Delta8,ht3_Delta9,ht3_r_EMIBD,ht3_F_Mean,ht3_F_SD,ht4_Wtheta,ht4_psi,ht4_psi2,ht4_Delta1,ht4_Delta2,ht4_Delta3,ht4_Delta4,ht4_Delta5,ht4_Delta6,ht4_Delta7,ht4_Delta8,ht4_Delta9,ht4_r_EMIBD,ht4_F_Mean,ht4_F_SD,ht5_Wtheta,ht5_psi,ht5_psi2,ht5_Delta1,ht5_Delta2,ht5_Delta3,ht5_Delta4,ht5_Delta5,ht5_Delta6,ht5_Delta7,ht5_Delta8,ht5_Delta9,ht5_r_EMIBD,ht5_F_Mean,ht5_F_SD" > r1_reference_table7.csv

for i in {1..5000}; do
  echo "Simulation $i"

  # Run SIMgen.py once and store output
  sim_output=$(python SIMgen.py)

  # Extract alpha and popsize from the single run
  transition_time=$(echo "$sim_output" | grep "Simulating with transition_time =" | awk '{print $NF}')
  alpha=$(echo "$sim_output" | grep "Simulating with alpha =" | awk '{print $NF}')
  popsize=$(echo "$sim_output" | grep "Simulating with popsize =" | awk '{print $NF}')

  # Rename generated VCFs (assuming they have generic names at first)
  mv EL2_SIM.vcf EL2_SIM_${i}.vcf
  mv JC69_SIM.vcf JC69_SIM_${i}.vcf
  mv TKF_SIM.vcf TKF_SIM_${i}.vcf
  mv JC69only_SIM.vcf JC69only_SIM_${i}.vcf

  #mv SIM.map SIM_${i}.map
  #mv SNP.map SNP_${i}.map

  # Compress/index
  bgzip -f EL2_SIM_${i}.vcf && bcftools index -f EL2_SIM_${i}.vcf.gz
  bgzip -f JC69_SIM_${i}.vcf && bcftools index -f JC69_SIM_${i}.vcf.gz
  bgzip -f TKF_SIM_${i}.vcf && bcftools index -f TKF_SIM_${i}.vcf.gz

  # Merge and sort
  bcftools concat EL2_SIM_${i}.vcf.gz JC69_SIM_${i}.vcf.gz TKF_SIM_${i}.vcf.gz -o merged_SIM_${i}.vcf
  bcftools sort merged_SIM_${i}.vcf -o merged_SIM_${i}.vcf

  # Generate distance matrix
  Rscript --no-save distancematrixgen.R merged_SIM_${i}.vcf $i

  # Run TDA & EMIBD9, and extract summary statistics
  export PATH="/data/proj2/home/students/c.yoon/Desktop/Thesis/EMIBD9/Bin:$PATH"
  find . -type f -name "EMIBD_${i}_*.par" | while read par_file; do
    echo "Running EM_IBD_P on $par_file"
    EM_IBD_P INP:$par_file
  done

  Rscript --no-save delta.R ${i} $i

  python tda.py "$transition_time" "$alpha" "$popsize" "$i"

  rm -f distanceMatrix_*
  rm -f seg*
  rm -f EMIBD_*
  rm -f merged_SIM_*
  rm -f Delta*
  rm -f F_*
  rm -f r_*


  # Generate distance matrix
  Rscript --no-save distancematrixgen.R JC69only_SIM_${i}.vcf $i

  # Run TDA & EMIBD9, and extract summary statistics
  export PATH="/data/proj2/home/students/c.yoon/Desktop/Thesis/EMIBD9/Bin:$PATH"
  find . -type f -name "EMIBD_${i}_*.par" | while read par_file; do
    echo "Running EM_IBD_P on $par_file"
    EM_IBD_P INP:$par_file
  done

  Rscript --no-save delta.R ${i} $i

  python tda.py "$transition_time" "$alpha" "$popsize" "$i"

  rm -f distanceMatrix_*
  rm -f Delta*
  rm -f F_*
  rm -f r_*
  rm -f seg*
  rm -f JC69*
  rm -f EL2*
  rm -f TKF*
  rm -f EMIBD_*
  rm -f merged_SIM_*
  rm -f hapibd_*
  rm -f SNP_*
done &> r1_reference_table7_output.log &


echo "Script is running in the background. Check 'r1_reference_table7_output.log' for progress."

##nohup bash r1_referencetable7.sh &> r1_referencetable7_output.log &
##pkill -f r1_referencetable7.sh
