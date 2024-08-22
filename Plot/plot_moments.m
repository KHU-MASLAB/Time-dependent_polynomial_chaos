function plot_moments(t_MC,t_PC,t_TDPC,mu_MC,mu_PC,mu_TDPC,tit)

close all;
fontsize=20;
titles = ["X disp.","Y disp","th angle","X vel.","Y vel.","th vel."];
endtime = t_MC(end);
h = (length(t_MC)-1)/endtime;
for i = 1:6
    minX = min(mu_MC(100:end,i));
    maxX = max(mu_MC(100:end,i));
    minmax = real([minX-(maxX-minX)*0.2,maxX+(maxX-minX)*0.2]);
    F = figure(i);
    hold on
    grid on
    plot(t_MC,mu_MC(:,i),'-k','LineWidth',1.5)
    plot(t_PC,mu_PC(:,i),'-.b','LineWidth',1.5)
    plot(t_TDPC,mu_TDPC(i,:),'--r','LineWidth',2)
    
    plot(t_MC,mu_MC(:,i),'squark','LineWidth',2,'MarkerIndices',1:h:length(t_MC),'MarkerSize',15)
    plot(t_PC,mu_PC(:,i),'ob','LineWidth',2,'MarkerIndices',1:h:length(t_PC),'MarkerSize',10)
    plot(t_TDPC,mu_TDPC(i,:),'^r','LineWidth',2,'MarkerIndices',1:h:length(t_TDPC),'MarkerSize',15)

    xlabel('time')
    ylabel(titles(i))
    legend('MC','PC','TD-PC')
    title(tit)
    set(gca,'FontName','Times New Roman','FontSize',fontsize);
%     set(F,'Position',[91 65 1500 850]);
    axis([t_MC(1),t_MC(end),minmax(1),minmax(2)])
end


