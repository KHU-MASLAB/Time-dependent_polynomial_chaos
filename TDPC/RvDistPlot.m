function RvDistPlot(rv)
dist = 'kernel';
close all
F = figure;
hold on 
grid on
for i = 1:size(rv,2)
    pd = fitdist(rv(:,i),dist,'Kernel','epanechnikov');
    Xsample = linspace(min(rv(:,i)),max(rv(:,i)),1000);
    Ysample = pdf(pd,Xsample);
    subplot(3,3,i)
    plot(Xsample,Ysample,'LineWidth',1.5)
end
set(gca,'FontName','Times New Roman','FontSize',20);
set(F,'Position',[91 65 1500 850]);