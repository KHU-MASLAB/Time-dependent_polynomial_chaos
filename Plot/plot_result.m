function plot_result(t_MC,t_PC,t_TDPC,mean_MC,mean_PC,mean_TDPC,var_MC,var_PC,var_TDPC,h)

%% mean of theta disp., theta vel.
close all
F1=figure;
hold on
grid on
plot(t_MC,mean_MC(:,3),'-k','LineWidth',1.5)
plot(t_PC,mean_PC(:,3),'-.b','LineWidth',1.5)
plot(t_TDPC,mean_TDPC(3,:),'--r','LineWidth',2)
plot(t_MC,mean_MC(:,3),'squark','LineWidth',2,'MarkerIndices',1:1/h:length(t_MC),'MarkerSize',15)
plot(t_PC,mean_PC(:,3),'ob','LineWidth',2,'MarkerIndices',1:1/h:length(t_PC),'MarkerSize',10)
plot(t_TDPC,mean_TDPC(3,:),'^r','LineWidth',2,'MarkerIndices',1:1/h:length(t_TDPC),'MarkerSize',15)
ax = gca;
ax.YTick = -2.5:0.25:-0.5;
legend('MC','PC','TD-PC')
xlabel('Time')
ylabel('Mean of Angular disp.')
set(gca,'FontName','Times New Roman','FontSize',20);
set(F1,'Position',[91 65 1500 850]);
axis([0,30,-2.5,-0.5])

F3 = figure;
grid on
hold on 
plot(t_MC,mean_MC(:,6),'-k','LineWidth',1.5)
plot(t_PC,mean_PC(:,6),'-.b','LineWidth',1.5)
plot(t_TDPC,mean_TDPC(6,:),'--r','LineWidth',2)
plot(t_MC,mean_MC(:,6),'squark','LineWidth',2,'MarkerIndices',1:1/h:length(t_MC),'MarkerSize',15)
plot(t_PC,mean_PC(:,6),'ob','LineWidth',2,'MarkerIndices',1:1/h:length(t_PC),'MarkerSize',10)
plot(t_TDPC,mean_TDPC(6,:),'^r','LineWidth',2,'MarkerIndices',1:1/h:length(t_TDPC),'MarkerSize',15)
legend('MC','PC','TD-PC')
xlabel('Time')
ylabel('Mean of Angular vel.')
set(gca,'FontName','Times New Roman','FontSize',20);
set(F3,'Position',[91 65 1500 850]);
axis([0,30,-2,2])
%% Variance of theta disp., theta vel.
F1=figure;
hold on
grid on
plot(t_MC,var_MC(:,3),'-k','LineWidth',1.5)
plot(t_PC,var_PC(:,3),'-.b','LineWidth',1.5)
plot(t_TDPC,var_TDPC(3,:),'--r','LineWidth',2)
plot(t_MC,var_MC(:,3),'squark','LineWidth',2,'MarkerIndices',1:1/h:length(t_MC),'MarkerSize',15)
plot(t_PC,var_PC(:,3),'ob','LineWidth',2,'MarkerIndices',1:1/h:length(t_PC),'MarkerSize',10)
plot(t_TDPC,var_TDPC(3,:),'^r','LineWidth',2,'MarkerIndices',1:1/h:length(t_TDPC),'MarkerSize',15)
ax = gca;
ax.YTick = 0 : 0.05 : 0.35;
legend('MC','PC','TD-PC')
xlabel('Time')
ylabel('Variance of Angular disp.')
set(gca,'FontName','Times New Roman','FontSize',20);
set(F1,'Position',[91 65 1500 850]);
axis([0,30,0,0.35])

F3 = figure;
grid on
hold on 
plot(t_MC,var_MC(:,6),'-k','LineWidth',1.5)
plot(t_PC,var_PC(:,6),'-.b','LineWidth',1.5)
plot(t_TDPC,var_TDPC(6,:),'--r','LineWidth',2)
plot(t_MC,var_MC(:,6),'squark','LineWidth',2,'MarkerIndices',1:1/h:length(t_MC),'MarkerSize',15)
plot(t_PC,var_PC(:,6),'ob','LineWidth',2,'MarkerIndices',1:1/h:length(t_PC),'MarkerSize',10)
plot(t_TDPC,var_TDPC(6,:),'^r','LineWidth',2,'MarkerIndices',1:1/h:length(t_TDPC),'MarkerSize',15)
legend('MC','PC','TD-PC')
xlabel('Time')
ylabel('Variance of Angular vel.')
set(gca,'FontName','Times New Roman','FontSize',20);
set(F3,'Position',[91 65 1500 850]);
ax = gca;
ax.YTick = 0:0.2:1.4;
axis([0,30,0,1.4])