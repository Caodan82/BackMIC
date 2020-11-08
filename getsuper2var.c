#include "mex.h"
#include "math.h"
void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int *vector_x,*best_c,*out1;
    int num_sample,avg,num1,num2,temc1;
    int i,k=1,if_fix=1;
    vector_x=mxGetData(prhs[0]);
    avg=*((int *)mxGetData(prhs[1]));
    num_sample=*((int *)mxGetData(prhs[2]));
    best_c=(int *) mxCalloc(num_sample,sizeof(int));
    for(i = 0;i < num_sample; i++){
        best_c[i]=0;
    }
    for(i=1;i<num_sample;i++){
        if (vector_x[i]!=vector_x[i-1]){
            if (if_fix==1){
                if ((i-best_c[k-1]) < avg){
                   num1= i-best_c[k-1];
                   temc1=i;
                   if_fix=0;
                }else{
                    best_c[k]=i;
                    k=k+1;
                    if_fix=1;
                }
            }else{
                num2=i-temc1;
                if ((num2+num1) < avg){
                    num1=num2+num1;
                    temc1=i;
                    if_fix=0;
                }else{
                    if ((avg-num1)>(num1+num2)-avg){
                        best_c[k]=i;
                        if_fix=1;
                        k=k+1;
                    }else{
                        best_c[k]=temc1;
                        k=k+1;
                        if (num2>=avg){
                            best_c[k]=i;
                            k=k+1;
                            if_fix=1;
                        }else{
                            num1=num2;
                            temc1=i;
                            if_fix=0;
                        }
                    }
                }
            }
         }
    }
    if (best_c[k-1]!=num_sample){
        best_c[k]=num_sample;
        plhs[0]=mxCreateNumericMatrix(1,k+1,mxINT32_CLASS,mxREAL);
        out1=mxGetData(plhs[0]);
        for (i=0;i<k+1;i++){
              out1[i]=best_c[i];
        }
    }else{
        plhs[0]=mxCreateNumericMatrix(1,k,mxINT32_CLASS,mxREAL);
        out1=mxGetData(plhs[0]);
        for (i=0;i<k;i++){
              out1[i]=best_c[i];
        }
    }
    mxFree(best_c);
    return;
}