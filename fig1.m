clear all

%read data
Pdata = xlsread('..\Figure4data.xls');

tdata = Pdata(:,1);
PCFdata = Pdata(:,2);
PCFstd = Pdata(:,3);

Adata = xlsread('..\Figure4data.xls',3);
Atdata= Adata(:,1);
ACFdata= Adata(:,2);
ACFstd= Adata(:,3);

figure(1)
hold on
errorbar(Atdata,ACFdata,ACFstd,'Ok','MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',1.5);
errorbar(tdata,PCFdata,PCFstd,'^g','MarkerSize',8,'MarkerEdgeColor','g','MarkerFaceColor','g','LineWidth',1.5);
plot(Atdata,0*Atdata,'k')
xlim([0,12])
ylim([-0.3,1])
xlabel('Generation (n)','interpreter','latex')
ylabel('Correlations','interpreter','latex')
legend('ACF exp','PCF exp')
set(gca,'FontSize',36)
pbaspect([1 1 1])