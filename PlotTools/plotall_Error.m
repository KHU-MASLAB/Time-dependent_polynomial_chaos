function plotall_Error(MC,PC,TDPC,tit)
t_MC = MC.t_MC; t_PC = PC.t_PC; t_TDPC = TDPC.t_TDPC;
eval("mu_MC = MC."+tit+"_MC;"); eval("mu_PC = PC."+tit+"_PC;"); eval("mu_TDPC = TDPC."+tit+"_TDPC;")

fontsize=20;
Error_PC = abs( (mu_PC - mu_MC(1:size(mu_PC,1),:))./mu_MC(1:size(mu_PC,1),:) );
Error_TDPC = abs( (mu_TDPC' - mu_MC(1:size(mu_TDPC',1),:))./mu_MC(1:size(mu_TDPC',1),:) );
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


