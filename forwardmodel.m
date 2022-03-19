function [Rh] = forwardmodel(x,i)
cd (strcat('./Rh_',num2str(i)));
save para.txt -ascii x
system('python ./calc_Rh_Peak.py');
Rh = importdata('Rh_sum.txt');
cd ('..')






