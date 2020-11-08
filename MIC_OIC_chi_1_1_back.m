function [MIC,thebestc1,thebestc2]=MIC_OIC_chi_1_1_back(data,B,c,n,threshold_a)
% OIC used the chi2test restrict segment by compare with two chi2test
% table (Local compare), the x_fix also use chi2test restrict (global compare)
% each table must campare with the 2*2 table; But all chi2test used adjust,
% include the df > 1
% n is the numbers of sample
% mutual_I_2 represent fixed the first column of data(x1)
% mutual_I_1 represent fixed the second column of data(x2)
% B=n^0.6
% c is max clumps (c*x)
%  see also MIC_2variable_chi.m : using Global chi2test compare
%           MIC_2variable_optimization2: using Local compare and Backtracking
%           MIC_2variable_optimization2_2: used the chi2test restrict segment by compare with two chi2test
%                table (Local compare),and replace log2(min(x,y)) with log2(x_fix)
%           MIC_OIC_chi_2: used the chi2test restrict segment by compare with two chi2test
%                        table (Local compare), Backtracking with Local
%                        compare
% MIC_OIC_chi_1_1_back: two  direction of  chi2test "MIC_OIC_chi_1_1"
% data=[Y,X];thebestc1,thebestc2分别是Y和X的剖分点
if nargin<5
    threshold_a=0.01;
end
mutual_I_2=zeros(round(B/2)-1,round(B/2)-1);
mutual_I_1=zeros(round(B/2)-1,round(B/2)-1);
mutual_chi_2=zeros(round(B/2)-1,round(B/2)-1);
mutual_chi_1=zeros(round(B/2)-1,round(B/2)-1);
last_c2=zeros(round(B/2)+1,round(B/2)+1);
last_c1=zeros(round(B/2)+1,round(B/2)+1);
% randnum=randperm(n)';
% data=data(randnum,:);
n=int32(n);
[value1,pos1]=sortrows(data,1);
[value2,pos2]=sortrows(data,2);
D1=nan(n,2);
D1=int32(D1);
for i=2:round(B/2)
    Q_x1=equipartitionYaxis2c(value1(:,1),int32(i),n);
    Q_x2=equipartitionYaxis2c(value2(:,2),int32(i),n);
    D1(pos1,1)=Q_x1;
    vector1=D1(pos2,1);
    avg=int32(n/(round(B/i)*c));
    c_x2=getsuper2var(vector1,avg,n);
    D1(pos2,2)=Q_x2;
    vector2=D1(pos1,2);
    c_x1=getsuper2var(vector2,avg,n);
    len1=int32(length(c_x1)-2);
    len2=int32(length(c_x2)-2);
    sub_max_seg=int32(round(B/i));
    chi2value=chi2inv(1-threshold_a,i-1);
    [temI2,tem_c2,chi2]=getmutualI2var_fix4(vector1,c_x2,sub_max_seg,int32(i),n,len2,chi2value);% equipartition x1
    [temI1,tem_c1,chi1]=getmutualI2var_fix4(vector2,c_x1,sub_max_seg,int32(i),n,len1,chi2value);% equipartition x2
    mutual_I_2(i-1,1:length(temI2))=temI2';
    mutual_chi_2(i-1,1:length(temI2))=chi2';
    last_c2(i-1,1:length(tem_c2))=tem_c2;
    mutual_I_1(1:length(temI1),i-1)=temI1;
    mutual_chi_1(1:length(temI1),i-1)=chi1;
    last_c1(1:length(tem_c1),i-1)=tem_c1';
end
% [ind_row,ind_col] = ind2sub(size(mutual_chi_2),find(mutual_chi_2~=0));
% pdist([ind_row,ind_col]);
last_c2_back=zeros(round(B/2)+1,round(B/2)+1);
last_c1_back=zeros(round(B/2)+1,round(B/2)+1);
mutual_I_2_back=zeros(round(B/2)-1,round(B/2)-1);
mutual_I_1_back=zeros(round(B/2)-1,round(B/2)-1);
mylastc1=cell(size(mutual_I_1,2),1);
mylastc2=cell(size(mutual_I_1,2),1);
sample_num=size(data,1);
for i=1:size(mutual_I_1,2)
    [~,pos]=max(mutual_I_1(:,i));
    bestc=last_c1(1:pos+1,i);
    if bestc(end)<sample_num
        bestc(length(bestc)+1)=sample_num;
    end
    mylastc1{i}=bestc(2:end-1);
    bestc=sort(bestc);
    for j=1:length(bestc)-1
        D1(pos1(bestc(j)+1:bestc(j+1)),1)=int32(j);
    end
    vector1=D1(pos2,1);
    avg=int32(n/(round(B/(length(bestc)-1))*c));
    c_x2=getsuper2var(vector1,avg,n);
    len2=int32(length(c_x2)-2);
    sub_max_seg=int32(round(B/(length(bestc)-1)));
    chi2value=chi2inv(1-threshold_a,length(bestc)-2);
    [temI2,tem_c2,~]=getmutualI2var_fix4(vector1,c_x2,sub_max_seg,int32(length(bestc)-1),n,len2,chi2value);% equipartition x1
    mutual_I_2_back(i,1:length(temI2))=temI2';
    last_c2_back(i,1:length(tem_c2))=tem_c2;
    [~,pos]=max(mutual_I_2(i,:));
    bestc=last_c2(i,1:pos+1);
    if bestc(end)<sample_num
        bestc(length(bestc)+1)=sample_num;
    end
    mylastc2{i}=bestc(2:end-1);
    bestc=sort(bestc);
    for j=1:length(bestc)-1
        D1(pos2(bestc(j)+1:bestc(j+1)),2)=int32(j);
    end
    vector2=D1(pos1,2);
    avg=int32(n/(round(B/(length(bestc)-1))*c));
    c_x1=getsuper2var(vector2,avg,n);
    len1=int32(length(c_x1)-2);
    sub_max_seg=int32(round(B/(length(bestc)-1)));
    chi2value=chi2inv(1-threshold_a,length(bestc)-2);
    [temI1,tem_c1,~]=getmutualI2var_fix4(vector2,c_x1,sub_max_seg,int32(length(bestc)-1),n,len1,chi2value);% equipartition x2
    mutual_I_1_back(1:length(temI1),i)=temI1;
    last_c1_back(1:length(tem_c1),i)=tem_c1';
end
MIC=max(max([mutual_I_1_back,mutual_I_2_back]));
position=find(mutual_I_2_back==MIC, 1);
if ~isempty(position)
    pos_row=mod(position(1)-1,round(B/2)-1)+1;
    pos_col=floor((position(1)-1)/(round(B/2)-1))+1;
    thebestc1=mylastc1{pos_row};% bestc of x1
    thebestc2=last_c2_back(pos_row,2:pos_col+1);% bestc of x2
    position2=find(mutual_I_1_back==MIC, 1);
    if ~isempty(position2)
        pos_row2=mod(position2(1)-1,round(B/2)-1)+1;
        pos_col2=floor((position2(1)-1)/(round(B/2)-1))+1;
        thebestc2_2=mylastc2{pos_col2};
        thebestc1_2=last_c1_back(2:pos_row2+1,pos_col2);
        if (length(thebestc2_2)+1)*(length(thebestc1_2)+1)<(length(thebestc2)+1)*(length(thebestc1)+1)
            thebestc1=thebestc1_2;
            thebestc2=thebestc2_2;
        end
    end
else
    position=find(mutual_I_1_back==MIC, 1);
    pos_row=mod(position(1)-1,round(B/2)-1)+1;
    pos_col=floor((position(1)-1)/(round(B/2)-1))+1;
    thebestc2=mylastc2{pos_col};
    thebestc1=last_c1_back(2:pos_row+1,pos_col);
end
end
