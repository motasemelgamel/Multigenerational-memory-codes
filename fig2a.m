clear all

%initialize
Nt=500;t=(0:1:Nt-1);PCFan=zeros(1,Nt);

%Noise strength
E=((0.2)^2)/2;Ip=((0.1)^2)/2;Is=((0.15)^2)/2;

%Parameters
beta=0.6;mu=0.3;

%Analytical ACF and PCF for size, sisters
syms x;
varan=2*(E+Is)/(beta*(2-beta)) + 2*Ip*((mu*(beta-1) -1)/(beta*(2-beta)*(mu^2 -1)*(1+mu*(beta-1))));
epsilon02=varan;y02=2*Ip/(1-mu^2);epsilon0y0=2*mu*Ip/((1-mu^2)*(1+mu*(beta-1)));
PCFan=limit(2*E*(1-(1-x).^(2*t))./(x*(2-x)) + epsilon02.*(1-x).^(2*t) + 2*epsilon0y0.*((1-x).^t).*(mu.^t - (1-x).^t)./(x+mu-1) + y02.*((((1-x).^t - mu.^t).^2)/(x+mu-1)^2),x,beta);
PCFannorm=PCFan/varan;
tauPan=1/(2*beta-beta^2) + 1/(2-2*mu) + 1/(2+2*mu) + 1/(1+mu*(beta-1)) + (Ip-Is*(mu^2 + mu -beta*mu -2))/(Ip*(mu*(beta-1)-1)+Is*(mu*(beta-1)+1)*(mu^2-1))

ACFan=limit(2*(E+Is)*((1-x).^t)/(x*(2-x)) - 2*Ip*(((1-x).^(t+1))/(x*(2-x)) + (mu.^(t+1))/(mu^2 -1))/((mu+x-1)*(1+mu*(x-1))),x,beta);
ACFnorm=ACFan/ACFan(1);
tauAan=(Ip*(1+mu*beta*(1-beta)+(beta-1)*mu^2)+E*(mu+1)*(1+mu*(beta-1))*(mu-1)^2 + Is*(mu+1)*(mu-1)^2*(1+mu*(beta-1)))/(beta*(mu-1)*(Ip*(mu*(beta-1)-1)+E*(mu*(beta-1)+1)*(mu^2-1)+Is*(mu*(beta-1)+1)*(mu^2-1)))

figure(1)
clf
hold on

plot(t,PCFannorm,'-^g','LineWidth',1.5)
xlim([0,12])
ylim([0,1])
plot(t,ACFnorm,'-Ok','LineWidth',1.5)
xlabel('Generation (n)','interpreter','latex')
ylabel('Correlations','interpreter','latex')

set(gca,'FontSize',36)
pbaspect([1 1 1])
legend('PCF','ACF')
title("Model 1 ACF and PCF for $\beta=$"+beta+", $\mu=$"+mu,'FontSize',36,'interpreter','latex')