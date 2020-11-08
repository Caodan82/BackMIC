#include "mex.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#define LEN sizeof(struct clumb)
struct clumb
{int subc;
 struct clumb *next;
};
struct clumb *insert(struct clumb *head,struct clumb *stud)
{
 struct clumb *sub_node0,*sub_node1,*sub_node2;  
 sub_node1=head;
 sub_node0=stud;
 sub_node2=0;
 if(head==0){
     head=sub_node0;
     sub_node0->next=0;
 }else{
     while((sub_node0->subc>sub_node1->subc)&&(sub_node1->next!=0)){
         sub_node2=sub_node1;
         sub_node1=sub_node1->next;         
     }
     if(sub_node0->subc<=sub_node1->subc){
         if(head==sub_node1) head=sub_node0;
         else sub_node2->next=sub_node0;
         sub_node0->next=sub_node1;   
     }else{
         sub_node1->next=sub_node0;
         sub_node0->next=0;    
     }
 }
 return(head);
}
struct clumb *del(struct clumb *head,int num)
{
 struct clumb *sub_node1,*sub_node2;
 sub_node1=head;
 if(head==0) mexPrintf("list is 0\n");
 else
 {while(num!=sub_node1->subc && sub_node1->next!=0){
     sub_node2=sub_node1;
     sub_node1=sub_node1->next;  
 }
 if(num==sub_node1->subc){
     if(sub_node1==head) head=sub_node1->next;
     else sub_node2->next=sub_node1->next;
     mxFree(sub_node1);
 }
 else mexPrintf("not been found\n"); 
}
 return(head);
}
void print(struct clumb *head)
{
 struct clumb *node;
 node=head;
 if(head!=0)
     do
     {
     mexPrintf("%d\n",node->subc);
     node=node->next;
     } while(node!=0);
}

void release(struct clumb *head)
{
 struct clumb *node;
 if(head!=0)
     do
     {
     node=head;
     head=head->next;
     mxFree(node);
     } while(head!=0);
}

void listtomatric(struct clumb *head,int *submat,int len_mat)
{
 struct clumb *node;
 int subcount=0;
 node=head;
 if(head!=0)
     do
     {
     if(subcount!=0) submat[subcount-1]=node->subc;
     subcount=subcount+1;
     node=node->next;
     } while(node!=0);
}
double myentropy(int mytable[],int tab_len,int num_sample)
{
    int i;
    double entr=0.0;
    for(i=0;i<tab_len;i++){
        if(mytable[i]!=0){
            entr=entr+(((double)mytable[i])/((double)num_sample))*(log(((double)mytable[i])/((double)num_sample))/log(2)); 
        }
    }
    entr=entr*(-1);
    return entr;
}
double mychi2test(int myT[], int myR[], int myC[], int lenT, int lenC, int lenR, int num_sample) {
    int i,j;
    double entr=0.0;
    for(i=0;i<lenR;i++){
        for(j=0;j<lenC;j++){
            entr=entr+pow(fabs((double)myT[j+i*lenC]-((double)myR[i]*(double)myC[j]/(double)num_sample))-0.5, 2)/((double)myR[i]*(double)myC[j]/(double)num_sample);
        }
    }
    return entr;
}
void mutual_I(double myI[],int vector_x[],struct clumb *head,int certain,int num_sample,int len_seg,int temcol)
{
    int i,j,subk=0,*subcertain,*temseg,*array_seg;
    int *array_seg_cer,*array_cer,*chi_array_seg_cer,*chi_array_cer,*chi_array_seg;
    double H1,H2,H3,mymulti;
    subcertain=(int *) mxCalloc(certain,sizeof(int));
    array_seg=(int *) mxCalloc(len_seg,sizeof(int));
    array_seg_cer=(int *) mxCalloc(certain*len_seg,sizeof(int));
    array_cer=(int *) mxCalloc(certain,sizeof(int));
    chi_array_seg=(int *) mxCalloc(2,sizeof(int));
    chi_array_seg_cer=(int *) mxCalloc(certain*2,sizeof(int));
    chi_array_cer=(int *) mxCalloc(certain,sizeof(int));
    for(i=0;i<certain*len_seg;i++){
        if(i<certain){
            array_cer[i]=0;
            chi_array_cer[i]=0;
            subcertain[i]=i+1;
        }
        if(i<len_seg) array_seg[i]=0;
        if(i<2) chi_array_seg[i]=0;
        if(i<certain*2) chi_array_seg_cer[i]=0;
        array_seg_cer[i]=0;
    }
    temseg=(int *) mxCalloc(len_seg,sizeof(int));
    listtomatric(head,temseg,len_seg);
    for(i=0;i<num_sample;i++){
        for(j=0;j<certain;j++){
            if(vector_x[i]==subcertain[j]){
                array_cer[j]=array_cer[j]+1;
                if(i<temseg[subk]){
                    array_seg_cer[subk+j*len_seg]=array_seg_cer[subk+j*len_seg]+1;
                    array_seg[subk]=array_seg[subk]+1;
                }else{
                    subk=subk+1;
                    array_seg[subk]=array_seg[subk]+1;
                    array_seg_cer[subk+j*len_seg]=array_seg_cer[subk+j*len_seg]+1;
                }
                if(len_seg>2){
                    if(subk>=temcol && subk<=(temcol+1)){
                        chi_array_cer[j]=chi_array_cer[j]+1;
                        if(subk==temcol){
                            chi_array_seg[0]=chi_array_seg[0]+1;
                            chi_array_seg_cer[0+j*2]=chi_array_seg_cer[0+j*2]+1;
                        }else{
                            chi_array_seg[1]=chi_array_seg[1]+1;
                            chi_array_seg_cer[1+j*2]=chi_array_seg_cer[1+j*2]+1;
                        }
                    }
                }
              break;
            }
        }
    }
    H1=myentropy(array_seg_cer,certain*len_seg,num_sample);
    H2=myentropy(array_seg,len_seg,num_sample);
    H3=myentropy(array_cer,certain,num_sample);
    if (len_seg<certain){
       myI[0]=(H2+H3-H1)/(log(len_seg)/log(2));
    }else{
       myI[0]=(H2+H3-H1)/(log(certain)/log(2));
    }
    if(len_seg>2){
            myI[1]=0.0;
            for(i=0;i<certain;i++){
                for(j=0;j<2;j++){
                    myI[1]=myI[1]+pow(fabs((double)chi_array_seg_cer[j+i*2]-((double)chi_array_cer[i]*(double)chi_array_seg[j]/((double)chi_array_seg[0]+(double)chi_array_seg[1])))-0.5,2)/((double)chi_array_cer[i]*(double)chi_array_seg[j]/((double)chi_array_seg[0]+(double)chi_array_seg[1]));
                }
            }
    }else{
       myI[1]=0.0; 
    }
    myI[2]=mychi2test(array_seg_cer,array_cer,array_seg,certain*len_seg,len_seg,certain,num_sample);
    mxFree(temseg);
    mxFree(subcertain);
    mxFree(array_cer);
    mxFree(array_seg);
    mxFree(array_seg_cer);
    mxFree(chi_array_cer);
    mxFree(chi_array_seg);
    mxFree(chi_array_seg_cer);
}
void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    int *vector_x,*best_c,*out2,mymark=0;
    int num_sample,segm,certain,len_bestc,len_segm,i,j,count,num,nofbestc,countc;
    double *myI,*chi2v,chi0,I0,*tem_I,mychivalue=0.0,*chivalue;
    struct clumb *head_mybestc,*head_tembestc,*node1,*node2,*node_tem;
    vector_x=mxGetData(prhs[0]);
    best_c=mxGetData(prhs[1]);
    segm=*((int *)mxGetData(prhs[2]));
    certain=*((int *)mxGetData(prhs[3]));
    num_sample=*((int *)mxGetData(prhs[4]));
    len_bestc=*((int *)mxGetData(prhs[5]));
    chivalue=mxGetPr(prhs[6]);
    if(len_bestc>segm-1) len_segm=segm-1;
    else len_segm=len_bestc;
    plhs[0]=mxCreateDoubleMatrix(len_segm,1,mxREAL);
    myI=mxGetPr(plhs[0]);
    tem_I=(double *) mxCalloc(3,sizeof(double));
    plhs[1]=mxCreateNumericMatrix(len_segm+2,1,mxINT32_CLASS,mxREAL);
    out2=mxGetData(plhs[1]);
    plhs[2]=mxCreateDoubleMatrix(len_segm,1,mxREAL);
    chi2v=mxGetPr(plhs[2]);
    head_mybestc=0;
    head_tembestc=0;
    node1=node2=(struct clumb *) mxMalloc(LEN);
    node1->subc=best_c[0];
    head_mybestc=node1;
    node2=node1;
    node1=(struct clumb *) mxMalloc(LEN);
    node1->subc=best_c[len_bestc+1];
    node2->next=node1;
    node1->next=0;  
    node1=node2=(struct clumb *) mxMalloc(LEN);
    for(i=1;i<len_bestc+1;i++){
        node1->subc=best_c[i];
        if(i==1) head_tembestc=node1;
        else node2->next=node1;
        node2=node1;
        node1=(struct clumb *) mxMalloc(LEN);
    }
    node2->next=0;    
    count=0;
    out2[0]=0;
    for(i=0;i<len_segm;i++){
        I0=0.0;
        chi0=0.0;
        node_tem=head_tembestc;
        if(head_tembestc!=0)
        do
        {
        j=node_tem->subc;
        node1=(struct clumb *) mxMalloc(LEN);
        node1->subc=j;
        head_mybestc=insert(head_mybestc,node1);
        countc=0;
        for(nofbestc=0;nofbestc<=i;nofbestc++){
            if(out2[nofbestc]<j){
                countc=countc+1;
            }           
        }
        countc=countc-1;
        tem_I[0]=0.0;
        tem_I[1]=0.0;
        tem_I[2]=0.0;
        mutual_I(tem_I,vector_x,head_mybestc,certain,num_sample,count+2,countc);
        if (fabs(tem_I[0])>fabs(I0)){
            I0=tem_I[0];
            num=j;
            mychivalue=tem_I[1];
            chi0=tem_I[2];
        }
        head_mybestc=del(head_mybestc,j);
        node_tem=node_tem->next;
        }while(node_tem!=0);
        if(i>0 && mychivalue<chivalue[0]){
            mymark=1;
            break;
        }
        myI[count]=I0;
        chi2v[count]=chi0;
        count++; 
        node1=(struct clumb *) mxMalloc(LEN);
        node1->subc=num;
        out2[i+1]=num;
        head_mybestc=insert(head_mybestc,node1);
        head_tembestc=del(head_tembestc,num);
    }
    if(mymark==1){
        out2[count+1]=num_sample;
    }else{
        out2[len_segm+1]=num_sample;
    }
    release(head_mybestc);
    release(head_tembestc);
    mxFree(tem_I);
    return;
}
