clear all

%initialize
Nt=500;t=(0:1:Nt-1);

%Noise strength
Ip=((0.1)^2)/2;Is=((0.3)^2)/2;

%Parameters
beta=0.5;lambda=-0.2;
r=(1+beta-lambda)/2;omega2=beta-r^2;

%Analytical ACF and PCF for size, sisters
am=(lambda-beta+1-sqrt((beta-lambda-1)^2 -4*lambda))/2;ap=(lambda-beta+1+sqrt((beta-lambda-1)^2 -4*lambda))/2;
varan=-2*(Ip*(1+ ap*am) + Is*(lambda^2 *(ap*am + 1)-2*lambda*(ap+am)+ am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
epsilon02=varan;delta02=-2*(2*Ip*(ap-1)*(am-1) + Is*beta^2 *(am*ap+1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));epsilon0delta0=2*(Ip*(ap-1)*(am-1) + Is*beta*(-lambda*(ap + am) + am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
PCFan=delta02*((am.^t - ap.^t)./(am - ap)).^2 + epsilon02*((am.^t *(am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap)).^2 + 2*epsilon0delta0*((am.^t - ap.^t)./(am-ap)).*((am.^t * (am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap));
PCFannorm=PCFan/varan;

ACFan=(2/((am-ap)*(am*ap -1))).*((  am.^(t)*(am*(Ip+Is)-(am^2 +1)*Is*lambda +am*Is*lambda^2)  )./(am^2 -1)  +  ( (ap).^(t)*(Is*lambda + ap^(2)*Is*lambda - ap*(Ip+Is+Is*lambda^2)) )./(ap^2 -1));
ACFnorm=ACFan/ACFan(1);


%plot results
figure(2)
clf
subplot(2,2,3)
plot(t,ACFnorm,'-ok','MarkerSize',8,'LineWidth',1.5)
hold on
plot(t,0*t,'k')
hold off
xlim([0,20])
ylim([-1,1])
title('ACF','interpreter','latex')
xlabel('Generation (n)','interpreter','latex')
set(gca,'FontSize',25)
pbaspect([1 1 1])
title("ACF for $r=$"+r+", $\omega^2=$"+omega2,'interpreter','latex')


%Parameters
beta=0.8;lambda=0.7;
r=(1+beta-lambda)/2;omega2=beta-r^2;

%Analytical ACF and PCF for size, sisters
am=(lambda-beta+1-sqrt((beta-lambda-1)^2 -4*lambda))/2;ap=(lambda-beta+1+sqrt((beta-lambda-1)^2 -4*lambda))/2;
varan=-2*(Ip*(1+ ap*am) + Is*(lambda^2 *(ap*am + 1)-2*lambda*(ap+am)+ am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
epsilon02=varan;delta02=-2*(2*Ip*(ap-1)*(am-1) + Is*beta^2 *(am*ap+1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));epsilon0delta0=2*(Ip*(ap-1)*(am-1) + Is*beta*(-lambda*(ap + am) + am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
PCFan=delta02*((am.^t - ap.^t)./(am - ap)).^2 + epsilon02*((am.^t *(am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap)).^2 + 2*epsilon0delta0*((am.^t - ap.^t)./(am-ap)).*((am.^t * (am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap));
PCFannorm=PCFan/varan;

ACFan=(2/((am-ap)*(am*ap -1))).*((  am.^(t)*(am*(Ip+Is)-(am^2 +1)*Is*lambda +am*Is*lambda^2)  )./(am^2 -1)  +  ( (ap).^(t)*(Is*lambda + ap^(2)*Is*lambda - ap*(Ip+Is+Is*lambda^2)) )./(ap^2 -1));
ACFnorm=ACFan/ACFan(1);


subplot(2,2,1)
plot(t,ACFnorm,'-ok','MarkerSize',8,'LineWidth',1.5)
hold on
plot(t,0*t,'k')
hold off
xlim([0,20])
ylim([-1,1])
title('ACF','interpreter','latex')
xlabel('Generation (n)','interpreter','latex')
set(gca,'FontSize',25)
pbaspect([1 1 1])
title("ACF for $r=$"+r+", $\omega^2=$"+omega2,'interpreter','latex')

%Parameters
beta=2;lambda=0.7;
r=(1+beta-lambda)/2;omega2=beta-r^2;

%Analytical ACF and PCF for size, sisters
am=(lambda-beta+1-sqrt((beta-lambda-1)^2 -4*lambda))/2;ap=(lambda-beta+1+sqrt((beta-lambda-1)^2 -4*lambda))/2;
varan=-2*(Ip*(1+ ap*am) + Is*(lambda^2 *(ap*am + 1)-2*lambda*(ap+am)+ am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
epsilon02=varan;delta02=-2*(2*Ip*(ap-1)*(am-1) + Is*beta^2 *(am*ap+1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));epsilon0delta0=2*(Ip*(ap-1)*(am-1) + Is*beta*(-lambda*(ap + am) + am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
PCFan=delta02*((am.^t - ap.^t)./(am - ap)).^2 + epsilon02*((am.^t *(am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap)).^2 + 2*epsilon0delta0*((am.^t - ap.^t)./(am-ap)).*((am.^t * (am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap));
PCFannorm=PCFan/varan;

ACFan=(2/((am-ap)*(am*ap -1))).*((  am.^(t)*(am*(Ip+Is)-(am^2 +1)*Is*lambda +am*Is*lambda^2)  )./(am^2 -1)  +  ( (ap).^(t)*(Is*lambda + ap^(2)*Is*lambda - ap*(Ip+Is+Is*lambda^2)) )./(ap^2 -1));
ACFnorm=ACFan/ACFan(1);


subplot(2,2,2)
plot(t,ACFnorm,'-ok','MarkerSize',8,'LineWidth',1.5)
hold on
plot(t,0*t,'k')
hold off
xlim([0,20])
ylim([-1,1])
title('ACF','interpreter','latex')
ylabel('ACF','interpreter','latex')
xlabel('Generation (n)','interpreter','latex')
set(gca,'FontSize',25)
pbaspect([1 1 1])
title("ACF for $r=$"+r+", $\omega^2=$"+omega2,'interpreter','latex')


%Parameters
beta=1;lambda=-0.4;
r=(1+beta-lambda)/2;omega2=beta-r^2;

%Analytical ACF and PCF for size, sisters
am=(lambda-beta+1-sqrt((beta-lambda-1)^2 -4*lambda))/2;ap=(lambda-beta+1+sqrt((beta-lambda-1)^2 -4*lambda))/2;
varan=-2*(Ip*(1+ ap*am) + Is*(lambda^2 *(ap*am + 1)-2*lambda*(ap+am)+ am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
epsilon02=varan;delta02=-2*(2*Ip*(ap-1)*(am-1) + Is*beta^2 *(am*ap+1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));epsilon0delta0=2*(Ip*(ap-1)*(am-1) + Is*beta*(-lambda*(ap + am) + am*ap +1))/((am^2 -1)*(ap^2 -1)*(am*ap -1));
PCFan=delta02*((am.^t - ap.^t)./(am - ap)).^2 + epsilon02*((am.^t *(am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap)).^2 + 2*epsilon0delta0*((am.^t - ap.^t)./(am-ap)).*((am.^t * (am+beta-lambda) - ap.^t *(ap+beta-lambda))./(am-ap));
PCFannorm=PCFan/varan;

ACFan=(2/((am-ap)*(am*ap -1))).*((  am.^(t)*(am*(Ip+Is)-(am^2 +1)*Is*lambda +am*Is*lambda^2)  )./(am^2 -1)  +  ( (ap).^(t)*(Is*lambda + ap^(2)*Is*lambda - ap*(Ip+Is+Is*lambda^2)) )./(ap^2 -1));
ACFnorm=ACFan/ACFan(1);


subplot(2,2,4)
plot(t,ACFnorm,'-ok','MarkerSize',8,'LineWidth',1.5)
hold on
plot(t,0*t,'k')
hold off
xlim([0,20])
ylim([-1,1])
title('ACF','interpreter','latex')
xlabel('Generation (n)','interpreter','latex')
set(gca,'FontSize',25)
pbaspect([1 1 1])
title("ACF for $r=$"+r+", $\omega^2=$"+omega2,'interpreter','latex')