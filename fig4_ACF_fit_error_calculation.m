clear all

%load experimental ACF and PCF values
Pdata = xlsread('..\Figure4data.xls');

tdata = Pdata(:,1);
PCFdata = Pdata(:,2);
PCFstd = Pdata(:,3);

Adata = xlsread('..\Figure4data.xls',3);
Atdata= Adata(:,1);
ACFdata= Adata(:,2);
ACFstd= Adata(:,3);

%read size and phase data
Nt=20;
[delta1,epsilon1]=Vashistha(Nt,'..\ALL_LB32_Size.xls','Expt1_1012');
[delta2,epsilon2]=Vashistha(Nt,'..\ALL_LB32_Size.xls','Expt1_1015');
[delta3,epsilon3]=Vashistha(Nt,'..\ALL_LB32_Size.xls','Expt3_2807');

delta=[delta3 delta2 delta1];
epsilon=[epsilon3 epsilon2 epsilon1];

k=1;
t_=(0:Nt-1);

Ip=((0.2)^2)/2;Is=((0.3)^2)/2;
am=@(beta,lambda) (lambda-beta+1-sqrt((beta-lambda-1).^2 -4*lambda))/2;ap=@(beta,lambda) (lambda-beta+1+sqrt((beta-lambda-1).^2 -4*lambda))/2;
ACFan=@(beta,lambda,t) (2./((am(beta,lambda)-ap(beta,lambda)).*(am(beta,lambda).*ap(beta,lambda) -1))).*((  am(beta,lambda).^(t).*(am(beta,lambda).*(Ip+Is)-(am(beta,lambda).^2 +1)*Is.*lambda +am(beta,lambda)*Is.*lambda.^2)  )./(am(beta,lambda).^2 -1)  +  ( (ap(beta,lambda)).^(t).*(Is*lambda + ap(beta,lambda).^(2)*Is.*lambda - ap(beta,lambda).*(Ip+Is+Is*lambda.^2)) )./(ap(beta,lambda).^2 -1));

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
        xs(k)=length(epsilon1);
        xs(k+1)=length(epsilon2);
        
        %ACF for individual lineages
     for n=1:length(epsilon1)
         A1(n)=mean((epsilon1(1:end-n+1)).*(epsilon1(n:end)))./mean((epsilon1(1:end)).*(epsilon1(1:end)));
         A2(n)=mean((epsilon2(1:end-n+1)).*(epsilon2(n:end)))./mean((epsilon2(1:end)).*(epsilon2(1:end)));
      end
        A=[A1' A2'];
        t1=(0:length(epsilon1)-1);
        fun1=@(x,t1) ACFan(x(2),x(1),t1)./ACFan(x(2),x(1),0);
        fun2=@(x,t1) ACFan(x(2),x(1),t1)./ACFan(x(2),x(1),0);
        fun=@(x,t1) [fun1(x,t1), fun2(x,t1)];
        options = optimoptions('lsqcurvefit','StepTolerance', 1e-16,'OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16);
        [x,resnorm,residual,exitflag,output,lambda,jacobian]=lsqcurvefit(fun,[0.01,0.01],t1',A,[],[],options);
        %error estimate
        residual=[residual(:,1) ; residual(:,2)];
        [Q,R] = qr(jacobian,0);
        mse = sum(abs(residual).^2)/(size(jacobian,1)-size(jacobian,2));
        Rinv = inv(R);
        Sigma = Rinv*Rinv'*mse;
        se(k,:) = sqrt(diag(Sigma));
        
        lambda_(k)=real(x(1));
        beta_(k)=real(x(2));
        k=k+1;
end

%transfer error to omega^2 and r space
se=full(se);

se_lambda=se(:,1);se_lambda=se_lambda(se_lambda~=0);
se_beta=se(:,2);se_beta=se_beta(se_beta~=0);

se_omega2=sqrt((0.5*(1-beta_'+lambda_')).^2 .*se_beta.^2 + (0.5*(1+beta_'-lambda_')).^2 .*se_lambda.^2);
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
