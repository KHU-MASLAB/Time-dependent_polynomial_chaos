function plot_convergance(con,tt,ii)

close all;
fontsize=20;
titles = ["X disp.","Y disp","$\theta$ angle","X vel.","Y vel.","$\theta$ angular vel."];
X = 1:6;
Y_PC = zeros(4,6);
Y_TDPC = zeros(4,6);
moment = ["mean","var","skew","kurt"];
for i = 1:6
    for j = 1:4
        eval("idt = find(con.t_PC"+num2str(i)+" == tt);");
        eval("Y_PC("+num2str(j)+","+num2str(i)+") = con.Error_"+moment(j)+"_PC"+num2str(i)+"(idt,ii);");
        eval("Y_TDPC("+num2str(j)+","+num2str(i)+") = con.Error_"+moment(j)+"_TDPC"+num2str(i)+"(idt,ii);");
    end    
end
for i = 1:4
    F=figure;
    semilogy(X,Y_PC(i,:),'-ob','LineWidth',1.5)
    hold on
    semilogy(X,Y_TDPC(i,:),'-^r','LineWidth',2)
    grid on
    xlabel('Order of basis')
    ylabel(moment(i))
    title(titles(ii))
    legend('PC','TD-PC')
    set(gca,'FontName','Times New Roman','FontSize',fontsize);
    set(F,'Position',[91 65 1500 850]);
    axis([1,6,10^-4,10^1])
    ax = gca;
    ax.XTick = 1:1:6;
end