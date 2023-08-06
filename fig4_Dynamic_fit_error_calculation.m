clear all

Nt=200;Ntpre=150;Nz=1000;t=(0:1:Nt-1);
Ip=((0.1)^2)/2;Is=((0.3)^2)/2;Nsample=200;Pall=zeros(Nsample,Nt);Aall=zeros(Nsample,Nt);

%read data
Nt=20;
[delta1,epsilon1]=Vashistha(Nt,'P:\Data\Vashistha 2021\ALL_LB32_Size.xls','Expt1_1012');
[delta2,epsilon2]=Vashistha(Nt,'P:\Data\Vashistha 2021\ALL_LB32_Size.xls','Expt1_1015');
[delta3,epsilon3]=Vashistha(Nt,'P:\Data\Vashistha 2021\ALL_LB32_Size.xls','Expt3_2807');

delta=[delta3 delta2 delta1];
epsilon=[epsilon3 epsilon2 epsilon1];
k=1;
for i=1:2:size(epsilon,2)-1
    clear epsilon1 epsilon2 delta1 delta2 A1 A2 A t1
    
    if length(epsilon{i}) == length(epsilon{i+1})
        epsilon1=(epsilon{i}-mean(epsilon{i}));
        epsilon2=(epsilon{i+1}-mean(epsilon{i+1}));
        elseif length(epsilon{i}) > length(epsilon{i+1})
        epsilon1=(epsilon{i}(1:length(epsilon{i+1}))-mean(epsilon{i}(1:length(epsilon{i+1}))));
        epsilon2=(epsilon{i+1}-mean(epsilon{i+1}));
        elseif length(epsilon{i}) < length(epsilon{i+1})
        epsilon1=(epsilon{i}-mean(epsilon{i}));
        epsilon2=(epsilon{i+1}(1:length(epsilon{i}))-mean(epsilon{i+1}(1:length(epsilon{i}))));
    end
    
    if length(delta{i}) == length(delta{i+1})
        delta1=(delta{i}-mean(delta{i}));
        delta2=(delta{i+1}-mean(delta{i+1}));
        elseif length(delta{i}) > length(delta{i+1})
        delta1=(delta{i}(1:length(delta{i+1}))-mean(delta{i}(1:length(delta{i+1}))));
        delta2=(delta{i+1}-mean(delta{i+1}));
        elseif length(delta{i}) < length(delta{i+1})
        delta1=(delta{i}-mean(delta{i}));
        delta2=(delta{i+1}(1:length(delta{i}))-mean(delta{i+1}(1:length(delta{i}))));
    end
        epsilon_=[epsilon1(1:end-1); epsilon2(1:end-1)];
        delta_=[delta1; delta2];
        X=[epsilon_(2:end),delta_(1:end-1)];
        Xresp=delta_(2:end);
        mdl=fitlm(X,Xresp,'Intercept',false);

        % estimate parameters beta and lambda based on linear fit of data
        betadest(k)=-mdl.Coefficients.Estimate(1);
        lambdadest(k)=mdl.Coefficients.Estimate(2);
        se_beta(k)=mdl.Coefficients.SE(1);
        se_lambda(k)=mdl.Coefficients.SE(2);
        k=k+1;
end

r_fit_dyn=(1+betadest-lambdadest)./2;omega2_fit_dyn=betadest-r_fit_dyn.^2;
se_omega2=sqrt((0.5*(1-betadest+lambdadest)).^2 .*se_beta.^2 + (0.5*(1+betadest-lambdadest)).^2 .*se_lambda.^2);
se_r=sqrt(0.25.*se_beta.^2 + 0.25.*se_lambda.^2);


function [delta,epsilon]=Vashistha(N,file,sheet_name)
%give the function file name with path
data=xlsread(file,sheet_name);

%extract lineage data into arrays
k=1;
for i=2:4:size(data,2)
    T(:,k)=data(:,i-1);
    T(:,k+1)=data(:,i-1);
    x(:,k)=data(:,i);
    x(:,k+1)=data(:,i+1);
    k=k+2;
end

%initialize birth size xn
xn(1,:)=x(1,:);
Tb(1,:)=T(1,:);

%Finding xn (birth size) of each generation and growth rate
for i=1:size(x,2)
    k=2;
    for j=3:size(x,1) 
        %Finding xn and Tn
        if ((x(j-1,i)>x(j-2,i) && x(j,i)<x(j-1,i)) && x(j,i)~=0 && x(j-1,i)~=0 && (abs(round(x(j-1,i)-x(j,i)))~=0))
            xn(k,i)=x(j,i);
            Tb(k,i)=T(j,i);Td(k-1,i)=T(j-1,i);
            xd(k-1,i)=x(j-1,i);
            k=k+1; %Update indecies
        end
    end
end

%Calculating phi_n and storing separate lineage data in cell arrays
for i=1:size(xn,2)
    xngen{:,i}=xn(xn(:,i)~=0,i);
    xdgen{:,i}=xd(xd(:,i)~=0,i);
    Tbgen{:,i}=Tb(Tb(2:end,i)~=0,i);
    Tdgen{:,i}=Td(Td(:,i)~=0,i);
    Tgen{:,i}=abs(Tdgen{:,i}-Tbgen{:,i});
    phigen{:,i}=log(xdgen{:,i}./xngen{:,i}(1:end-1,:));
    alpha_{:,i}=phigen{:,i}./Tgen{:,i};
end

%N is Lineage Length
k=1;
for j=1:2:size(xngen,2)-1
    if (length(xngen{j})>=N && length(xngen{j+1})>=N)
    epsilon{k}=log(xngen{j}(1:end)/mean(xngen{j}));
    epsilon{k+1}=log(xngen{j+1}(1:end)/mean(xngen{j+1}));
    delta{k}=(phigen{j}-mean(phigen{j}));
    delta{k+1}=(phigen{j+1}-mean(phigen{j+1}));
    phi{k}=phigen{j};
    phi{k+1}=phigen{j+1};
    k=k+2;
    end
end
end

b=((-1+beta-lambda).^2 -4*lambda)./4;
c=1/2 - beta./2 + lambda./2;
a1=c -sqrt(b);
a2=c +sqrt(b);

if((abs(a1)<1  & b<0) | (abs(a1)<1 & abs(a2)<1 & b>0))
a=1;
else
a=0;
end
end