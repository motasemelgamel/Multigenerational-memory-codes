clear all

%Noise strength
E=((0.2)^2)/2;Ip=((0.1)^2)/2;Is=((0.15)^2)/2; %(standard deviation)^2 /2

figure(1)
clf
hold on

k1=1;
for mu=0.05:0.01:0.95
k2=1;
for beta=0.05:0.01:1.95
tauPan=1/(2*beta-beta^2) + 1/(2-2*mu) + 1/(2+2*mu) + 1/(1+mu*(beta-1)) + (Ip-Is*(mu^2 + mu -beta*mu -2))/(Ip*(mu*(beta-1)-1)+Is*(mu*(beta-1)+1)*(mu^2-1));
tauAan=(Ip*(1+mu*beta*(1-beta)+(beta-1)*mu^2)+E*(mu+1)*(1+mu*(beta-1))*(mu-1)^2 + Is*(mu+1)*(mu-1)^2*(1+mu*(beta-1)))/(beta*(mu-1)*(Ip*(mu*(beta-1)-1)+E*(mu*(beta-1)+1)*(mu^2-1)+Is*(mu*(beta-1)+1)*(mu^2-1)));
PInfty(k1,k2)=E/(E+Is+Ip*((mu*(beta-1)-1)/((mu^2 -1)*(mu*(beta-1)+1))));
DeltaTau(k1,k2)=tauPan-tauAan;
k2=k2+1;
end
scatter(PInfty(k1,:),DeltaTau(k1,:),(3*mu)^3,(0.05:0.01:1.95),'filled')
k1=k1+1;
end

plot(0,2,'kp','MarkerSize',30)
text(0.03,2,'Experimental Observation','VerticalAlignment','bottom','HorizontalAlignment','left')
colormap(jet)
bar=colorbar;
title("Difference between PCF and ACF fall off time scale plotted against PCF asymptote. $I_e$ = "+E+", $I_p$= "+Ip+", $\frac{I_e}{I_p}$= "+(E/Ip),'interpreter','latex')
xlabel('$P_\infty$','interpreter','latex')
ylabel('$\Delta \tau$','interpreter','latex')
xlim([0 0.7])
ylim([-10 10])
set(gca,'FontSize',36)
pbaspect([1 1 1])
yline(0)
L(1) = plot(nan, nan, '.k','MarkerSize',(3*0.95)^3);L(2) = plot(nan, nan, '.k','MarkerSize',(3*0.7)^3);L(3) = plot(nan, nan, '.k','MarkerSize',(3*0.55)^3);
legend(L, {'$\mu=0.95$','$\mu=0.7$','$\mu=0.55$'},'interpreter','latex');
set(get(bar,'Title') ,'String','$\beta$','interpreter','latex');