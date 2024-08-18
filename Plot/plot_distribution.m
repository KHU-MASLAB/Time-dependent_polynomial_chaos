function plot_distribution(t_PC,t_TDPC,thSample_MC,thSample_PC,thSample_TDPC,tt,variable,filename,Distribution)
dist = 'kernel';
idt1 = tt+1;
Np_TDPC = size(thSample_TDPC,1);

pd1 = fitdist(datasample(thSample_MC(:,idt1),Np_TDPC),dist,'Kernel','epanechnikov');
pd2 = fitdist(thSample_PC(:,idt1),dist,'Kernel','epanechnikov');
pd3 = fitdist(thSample_TDPC(:,idt1),dist,'Kernel','epanechnikov');
Xsample = linspace(min(thSample_MC(:,idt1)),max(thSample_MC(:,idt1)),1000);
Ysample1 = pdf(pd1,Xsample);
Ysample2 = pdf(pd2,Xsample);
Ysample3 = pdf(pd3,Xsample);
F = figure;
hold on 
grid on
plot(Xsample,Ysample1,'k','LineWidth',1.5)
plot(Xsample,Ysample2,'-.b','LineWidth',1.5)
plot(Xsample,Ysample3,'--r','LineWidth',2)
plot(Xsample,Ysample1,'squarek','LineWidth',2,'MarkerIndices',1:50:1000,'MarkerSize',15)
plot(Xsample,Ysample2,'ob','LineWidth',2,'MarkerIndices',1:50:1000,'MarkerSize',10)
plot(Xsample,Ysample3,'^r','LineWidth',2,'MarkerIndices',1:50:1000,'MarkerSize',15)
% histogram(datasample(thSample_MC(:,idt1),Np_TDPC),80,'Normalization','pdf','FaceColor','k')
% histogram(thSample_PC(:,idt2),80,'Normalization','pdf','FaceColor','b')
% histogram(thSample_TDPC(:,idt3(1)),80,'Normalization','pdf','FaceColor','r')
axis([min(thSample_MC(:,idt1)),max(thSample_MC(:,idt1)),0,max(Ysample1)*1.2])

legend('MC','PC','TD-PC')
xlabel(variable)
ylabel("f(\xi)")
set(gca,'FontName','Times New Roman','FontSize',20);
set(F,'Position',[91 65 1500 850]);
saveas(F,"plot/"+Distribution+"/PDF_time" + tt + "_"+variable+"_"+filename +".png")
% ax = gca;
% ax.YTick = 0:0.2:1.4;
% axis([0,20,0,1.4])
% histogram(th_sample_MC(:,idt1),50,'FaceColor','k','Normalization','pdf')