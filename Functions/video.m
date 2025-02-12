function video(t,disp,Prop)

lw = 2;
close all
F1 = figure();
hold on; grid on;
% set(F1,'Position',[90 65 1500 850]); 
axis([-2,2,-2,2])
pbaspect([1 1 1])
n = 1;
for i = 1:(1/Prop.h)/10:length(t)
    cla;
    plot([0,disp(i,1)*2],[0,disp(i,2)*2],'-ok','linewidth', lw)
    
    str = "time : " + num2str(t(i));
    txt = text(-1,1.2,str);
    txt.FontSize = 25;
    K2(n) = getframe;
    n = n+1;
    pause(0.1)
end

filename = 'SinglePendulum.gif';
for n=1:length(K2)
    im = frame2im(K2(n));
    [imind,cm] = rgb2ind(im,256);
      if n == 1
          imwrite(imind,cm,filename,'gif','DelayTime',0, 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif','DelayTime',0,'WriteMode','append');
      end
end