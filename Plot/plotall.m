function plotall(t_MC,t_PC,t_TDPC,mean_MC,mean_PC,mean_TDPC,tit)

close all;
fontsize=20;
titles = ["X disp.","Y disp","th angle","X vel.","Y vel.","th vel."];
endtime = t_MC(end);
h = (length(t_MC)-1)/endtime;
for i = 1:6
    minX = min(mean_MC(100:end,i));
    maxX = max(mean_MC(100:end,i));
    minmax = real([minX-(maxX-minX)*0.2,maxX+(maxX-minX)*0.2]);
    F=figure;
    hold on
    grid on
    plot(t_MC,mean_MC(:,i),'-k','LineWidth',1.5)
    plot(t_PC,mean_PC(:,i),'-.b','LineWidth',1.5)
    plot(t_TDPC,mean_TDPC(i,:),'--r','LineWidth',2)
    
    plot(t_MC,mean_MC(:,i),'squark','LineWidth',2,'MarkerIndices',1:h:length(t_MC),'MarkerSize',15)
    plot(t_PC,mean_PC(:,i),'ob','LineWidth',2,'MarkerIndices',1:h:length(t_PC),'MarkerSize',10)
    plot(t_TDPC,mean_TDPC(i,:),'^r','LineWidth',2,'MarkerIndices',1:h:length(t_TDPC),'MarkerSize',15)

    xlabel('time')
    ylabel(titles(i))
    legend('MC','PC','TD-PC')
    title(tit)
    set(gca,'FontName','Times New Roman','FontSize',fontsize);
    set(F,'Position',[91 65 1500 850]);
    axis([t_MC(1),t_MC(end),minmax(1),minmax(2)])
end


