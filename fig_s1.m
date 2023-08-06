clear all

lambda1= 0.9.*rand(1500,1);
beta1=(2*lambda1 + 1).*rand(1500,1);
Nz=1;
Nt=20;t=(0:1:Nt-1);Ntpre=500;epsilon=zeros(Ntpre,Nz);epsilon1=zeros(Nt+Ntpre,Nz);epsilon2=zeros(Nt+Ntpre,Nz);
Ip=((0.2)^2);Is=((0.3)^2);delta=zeros(Ntpre,Nz);delta1=zeros(Nt+Ntpre,Nz);delta2=zeros(Nt+Ntpre,Nz);

am=@(beta,lambda) (lambda-beta+1-sqrt((beta-lambda-1).^2 -4*lambda))/2;ap=@(beta,lambda) (lambda-beta+1+sqrt((beta-lambda-1).^2 -4*lambda))/2;
ACFan=@(beta,lambda,t) (2./((am(beta,lambda)-ap(beta,lambda)).*(am(beta,lambda).*ap(beta,lambda) -1))).*((  am(beta,lambda).^(t).*(am(beta,lambda).*(Ip+Is)-(am(beta,lambda).^2 +1)*Is.*lambda +am(beta,lambda)*Is.*lambda.^2)  )./(am(beta,lambda).^2 -1)  +  ( (ap(beta,lambda)).^(t).*(Is*lambda + ap(beta,lambda).^(2)*Is.*lambda - ap(beta,lambda).*(Ip+Is+Is*lambda.^2)) )./(ap(beta,lambda).^2 -1));
ACFan1=@(beta,lambda) (2./((am(beta,lambda)-ap(beta,lambda)).*(am(beta,lambda).*ap(beta,lambda) -1))).*((  am(beta,lambda).^(0).*(am(beta,lambda).*(Ip+Is)-(am(beta,lambda).^2 +1)*Is.*lambda +am(beta,lambda)*Is.*lambda.^2)  )./(am(beta,lambda).^2 -1)  +  ( (ap(beta,lambda)).^(0).*(Is*lambda + ap(beta,lambda).^(2)*Is.*lambda - ap(beta,lambda).*(Ip+Is+Is*lambda.^2)) )./(ap(beta,lambda).^2 -1));

r_=(1+beta1-lambda1)./2;omega_=real(sqrt(beta1-r_.^2));omega=omega_(omega_~=0);r=r_(omega_~=0);lambda=lambda1(omega_~=0);beta=beta1(omega_~=0);

%Generate data using model
for k=1:length(omega)
    lambda_=lambda(k);
    beta_=beta(k);
    for q=1:Ntpre
        eta=sqrt(Is)*randn;zeta=sqrt(Ip)*randn; %draw noises
        epsilon(q+1)=epsilon(q)+delta(q)+eta;
        delta(q+1)=-beta_*epsilon(q+1)+zeta;
    end
        epsilon1(1:Ntpre+1)=epsilon;delta1(1:Ntpre+1)=delta;
        epsilon2(1:Ntpre+1)=epsilon;delta2(1:Ntpre+1)=delta;
    for p=1:Nt
        eta1=sqrt(Is)*randn;zeta1=sqrt(Ip)*randn;eta2=sqrt(Is)*randn;zeta2=sqrt(Ip)*randn;
        epsilon1(Ntpre+p+1)=epsilon1(Ntpre+p)+delta1(Ntpre+p)+eta1;
        delta1(Ntpre+p+1)=-beta_*epsilon1(Ntpre+p+1)+lambda_*delta1(Ntpre+p)+zeta1;
        epsilon2(Ntpre+p+1)=epsilon2(Ntpre+p)+delta2(Ntpre+p)+eta2;
        delta2(Ntpre+p+1)=-beta_*epsilon2(Ntpre+p+1)+lambda_*delta2(Ntpre+p)+zeta2;
    end
epsilon1_=epsilon1(Ntpre+1:Ntpre+1+Nt);
epsilon2_=epsilon2(Ntpre+1:Ntpre+1+Nt);
for n=1:Nt
    %Size ACF, sisters
    A1(n)= mean((epsilon1_(1:end-n+1)).*(epsilon1_(n:end)))./mean((epsilon1_(1:end)).*(epsilon1_(1:end)));
    A2(n)= mean((epsilon2_(1:end-n+1)).*(epsilon2_(n:end)))./mean((epsilon2_(1:end)).*(epsilon2_(1:end)));
end

	A=[A1' A2'];
    fun1=@(x,t1) ACFan(x(2),x(1),t1)./ACFan1(x(2),x(1));
    fun=@(x,t1) [fun1(x,t1), fun1(x,t1)];
    options = optimoptions('lsqcurvefit','StepTolerance', 1e-16,'OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16);
    x=lsqcurvefit(fun,[0.01,0.01],t',A,[],[],options);
    lambda_fit(k)=real(x(1));
    beta_fit(k)=real(x(2));
end

r_fit=(1+beta_fit-lambda_fit)./2;omega_fit=real(sqrt(beta_fit-r_fit.^2));

percentage=sum(omega_fit&omega')/length(omega) *100;

figure(1)
clf
hold on
plot(omega,omega_fit,'O')
plot(omega,omega,'k')
xlabel('Actual frequency, \omega (gen.^{-1})')
ylabel('Fitted frequency, \omega (gen.^{-1})')
title("N="+Nt+" gen., "+"off-axes="+percentage+"%")
set(gca,'FontSize',30)
pbaspect([1 1 1])
