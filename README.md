# BackMIC
# A new algorithm for the calculation of the Maximal Information Coeffecient(MIC)
# BEFORE STARTING
# Make sure that c++ has installed in your computer;
# Make sure Matlab has installed in your computer.
# 1.  Down load all the files.equipartitionYaxis2c.c getmutualI2var_fix4.c getsuper2var.c make_mex.m MIC_back.m 
# 2.  run“make_mex.m” to compile the equipartitionYaxis2c.c, getsuper2var.c and getmutualI2var_fix4.c to mex files;
# 3.  Preparing a paired variable X and Y  to measure their correlation.   X and Y are both column vectors with the same sample size n. 
# 4.  Before computing the MIC of X and Y with BackMIC, you need first scramble the samples, the codes are as follows:
#    data=[X,Y];
#    num=randperm(size(data,1));
#    data=data(num',:); 
# 5. Run "MIC_back.m”to obtain the MIC value of X and Y.
# [BackMIC,thebestc1,thebestc2]=MIC_back(data,B(n),c,n,threshold ); 
# n is the sample size; 
# B(n) is always set to  n0.6;
# c is always set to 5; 
# threshold is always set to 0.01.
# thebestc1 and thebestc2 are the segment point along the x-axis and y-axis.
