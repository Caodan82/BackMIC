#include "mex.h"
#include "math.h"
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    int segment, num_sample, avg, num1, num2, temc1, temsegment;
    int i, j, k=1, if_fix=1;
    int *out1, *best_c;
    double *vector_x;
    vector_x=mxGetPr(prhs[0]);
    segment=*((int *)mxGetData(prhs[1]));
    num_sample=*((int *)mxGetData(prhs[2]));
    plhs[0]=mxCreateNumericMatrix(num_sample, 1, mxINT32_CLASS, mxREAL);
    out1=mxGetData(plhs[0]);
    best_c=(int *) mxCalloc(segment+1, sizeof(int));
    for(i = 0;i <= segment; i++){
        best_c[i]=0;
    }
    temsegment=1;
    avg=num_sample/segment;
    for(i=1;i<num_sample;i++){
        if ((num_sample-best_c[k-1])<=(segment-temsegment+1)){
            segment=temsegment;
        }
        if (temsegment==segment){
            if (if_fix==1){
                out1[i-1]=segment;
            }else{
                for(j=best_c[k-1];j<i;j++){
                    out1[j]=temsegment;
                }
                if_fix=1;
            }
            if (i==num_sample-1){
                out1[i]=temsegment;
            }
        }else{
            if (vector_x[i]!=vector_x[i-1]){
                if (if_fix==1){
                    if ((i-best_c[k-1]) < avg){
                        num1= i-best_c[k-1];
                        temc1=i;
                        if_fix=0;
                    }else{
                        best_c[k]=i;
                        if_fix=1;
                        for(j=best_c[k-1];j<i;j++){
                            out1[j]=temsegment;
                        }
                        avg=(num_sample-i)/(segment-temsegment);
                        temsegment=temsegment+1;
                        k=k+1;
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
                            for(j=best_c[k-1];j<i;j++){
                                out1[j]=temsegment;
                            }
                            avg=(num_sample-i)/(segment-temsegment);
                            temsegment=temsegment+1;
                            k=k+1;
                            num1=0;
                        }else{
                            best_c[k]=temc1;
                            for(j=best_c[k-1];j<temc1;j++){
                                out1[j]=temsegment;
                            }
                            avg=(num_sample-temc1)/(segment-temsegment);
                            temsegment=temsegment+1;
                            k=k+1;
                            num1=num2;
                            temc1=i;
                            if_fix=0;
                        }
                    }
                }
            }else{
                if (if_fix==1 && i==num_sample-1){
                    for(j=best_c[k-1];j<num_sample;j++){
                        out1[j]=temsegment;
                    }
                } 
            }
        }
    }
    if (if_fix==0){
        if (num1==0){
            for (j=best_c[k-1];j<num_sample;j++){
                out1[j]=temsegment;
            }
        }else{
            num2=num_sample-temc1;
            if ((num2+num1) < avg){
                for (j=best_c[k-1];j<num_sample;j++){
                    out1[j]=temsegment;
                }
            }else{
                if ((avg-num1)>(num1+num2)-avg){
                    for (j=best_c[k-1];j<num_sample;j++){
                        out1[j]=temsegment;
                    }
                }else{
                    for (j=best_c[k-1];j<best_c[k-1]+num1;j++){
                        out1[j]=temsegment;
                    }
                    for (j=best_c[k-1]+num1;j<num_sample;j++){
                        out1[j]=temsegment+1;
                    }
                }
            }
        }
    }
    mxFree(best_c);
    return;
}
