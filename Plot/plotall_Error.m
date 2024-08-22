function plotall_Error(t_PC,t_TDPC,mean_MC,mean_PC,mean_TDPC,tit)

fontsize=20;
Error_PC = abs( (mean_PC - mean_MC(1:size(mean_PC,1),:))./mean_MC(1:size(mean_PC,1),:) );
Error_TDPC = abs( (mean_TDPC' - mean_MC(1:size(mean_TDPC',1),:))./mean_MC(1:size(mean_TDPC',1),:) );
titles = ["X disp.","Y disp","theta angle","X vel.","Y vel.","theta angular vel."];
endtime = t_TDPC(end);
h = (length(t_TDPC)-1)/endtime;
for i = 1:6
    F = figure(i+6);
    semilogy(t_PC,Error_PC(:,i),'-.b','LineWidth',1.5)
    hold on
    semilogy(t_TDPC,Error_TDPC(:,i),'--r','LineWidth',2)
    semilogy(t_PC,Error_PC(:,i),'ob','LineWidth',2,'MarkerIndices',1:h:length(t_PC),'MarkerSize',10)
    semilogy(t_TDPC,Error_TDPC(:,i),'^r','LineWidth',1.5,'MarkerIndices',1:h:length(t_TDPC),'MarkerSize',10)    
    grid on
    xlabel('time')
    ylabel(titles(i))
    legend('PC','TD-PC')
    title(tit)
    set(gca,'FontName','Times New Roman','FontSize',fontsize);
%     set(F,'Position',[91 65 1500 850]);
    axis([0,30,10E-6,10E2])
end


