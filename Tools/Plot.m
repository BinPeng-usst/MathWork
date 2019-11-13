clear all; clc; close all;
%% 3D plot
% fv0=[0.25
% 0.26
% 0.34
% 0.34
% 0.36
% 0.36
% ];
% s=[0.42
% 0.43
% 0.54
% 0.55
% 0.52
% 0.54
% ];
% tested=[0.34
% 0.74
% 0.86
% 0.66
% 0.72
% 0.67
% ];
% GP1=[0.63
% 0.5
% 0.68
% 0.75
% 0.72
% 0.77
% ];
% plot3(fv0, s,tested,'-.','linewidth',2);
% hold on; 
% plot3(fv0, s, GP1, '-','linewidth',2);
% grid on;
% box off;

%% PDF及概率面积图
Tested=0.74;
mu=0.49;
sigma=0.25;
pd = makedist('Normal',mu,sigma);
XStart=0;
XEnd=2.0;
Ylim=7;

dX=linspace(XStart,XEnd,1000);
P_dX=linspace(XStart,Tested,1000);


Ar=cdf(pd,Tested)*100;
SAr=num2str(Ar,'%.0f%%');

Smu=num2str(mu,'%.2f');
Ssigma=num2str(sigma,'%.2f');


plot(dX,pdf(pd,dX),'r','LineWidth',2);
hold on;
area(P_dX, pdf(pd,P_dX));
hold on;
plot(Tested, pdf(pd,Tested), 'r*','MarkerSize',8);
hold on;
% annotation('textarrow',[(Tested+0.5)/(XEnd-XStart) Tested/(XEnd-XStart)],[(pdf(pd,Tested)+0.5)/(XEnd-XStart) pdf(pd,Tested)/Ylim],'String',SAc,'FontName','Timesnewroman','FontSize',17)
 text(Tested, pdf(pd,Tested),['\leftarrow'],'Rotation',135,'FontName','Timesnewroman','FontSize',23); 
 hold on;
 str={'experimental', 'result'};
 text(Tested-0.3, pdf(pd,Tested)+1.5,str,'FontName','Timesnewroman','FontSize',23); 
 hold on
 str={['\it\mu=','\rm',Smu,' ,\it\sigma=','\rm',Ssigma],['Assuring rate->',SAr]};
 text(0.85, 4.5,str,'FontName','Timesnewroman','FontSize',23); 
 hold on
set(gca,'FontName','Timesnewroman','FontSize',23,'XAxisLocation','origin','xlim',[XStart XEnd],'xtick',[XStart:0.5: XEnd],'ylim',[0,Ylim],'ytick',[0:1.0:Ylim]);
box off;
 ylabel('Probability Density of \it\tau ','FontName','Timesnewroman','FontSize',23);
 xlabel('\it\tau \rm (N/mm^2)','FontName','Timesnewroman','FontSize',23);
%  print('-f1','H:\MathWork\Figures\v6','-painters','-dmeta','-r600');


%% 方差均值图
% No=[1:1:6];
% R=[0.74
% 0.86
% 0.66
% 0.67
% 0.72
% 0.34
% ];
% Mean=[0.5398
% 0.7003
% 0.7007
% 0.712
% 0.7034
% 0.51
% ];
% Sigma=[0.1296
% 0.1223
% 0.1232
% 0.1225
% 0.1222
% 0.1295
% ];
% errorbar(No,Mean,Sigma);
% hold on;
% plot(No,R);
%  set(gca,'FontName','Timesnewroman','FontSize',35, 'xlim',[0,7],'xtick',[0:1:7],'ylim',[0,1],'ytick',[0:0.2:1]);

%% subplot(1,2,2);%玫瑰图
% ypred=[0.10 
% 0.20 
% 0.18 
% 0.13 
% 0.38 
% 0.24 
% 0.67 
% 0.40 
% 0.23 
% 0.25 
% 0.23 
% 0.46 
% 0.56 
% 0.51 
% 0.63 
% 0.63 
% 0.65 
% 0.64 
% 0.49 
% ];
% Gama12=[0.12 
% 0.18 
% 0.17 
% 0.16 
% 0.26 
% 0.23 
% 0.65 
% 0.39 
% 0.21 
% 0.23 
% 0.20 
% 0.27 
% 0.27 
% 0.44 
% 0.56 
% 0.57 
% 0.59 
% 0.58 
% 0.42 
% ];
% Gama135=[0.13 
% 0.19 
% 0.17 
% 0.16 
% 0.27 
% 0.24 
% 0.67 
% 0.40 
% 0.22 
% 0.23 
% 0.21 
% 0.28 
% 0.28 
% 0.45 
% 0.58 
% 0.58 
% 0.60 
% 0.59 
% 0.43 
% ];
% CollectO2=[0.09 
% 0.21 
% 0.16 
% 0.11 
% 0.35 
% 0.28 
% 0.76 
% 0.41 
% 0.31 
% 0.42 
% 0.30 
% 0.66 
% 0.68 
% ];
% TestO=[0.74 
% 0.86 
% 0.66 
% 0.67 
% 0.72 
% 0.34 
% ];
% 
% s1=size(ypred,1);
% s2=size(CollectO2,1);
% s3=size(TestO,1);
% delta=2*pi/s1;
% polarplot([0:delta:(2*pi-delta)],ypred,'r-^','linewidth',3,'markersize',8);
% hold on;
% % polarplot([0:delta:(2*pi-delta)],Gama12,'b:','linewidth',2,'markersize',8);
% % hold on;
% polarplot([0:delta:(2*pi-delta)],Gama135,'k-.','linewidth',2,'markersize',8);
% hold on;
% polarscatter([0:delta:(s2-1)*delta]',CollectO2,300,'b');
% hold on;
% polarscatter([s2*delta:delta:(2*pi-delta)]',TestO,300,'p','r','filled');
% hold on;
% TTick=[0:rad2deg(delta):rad2deg(2*pi-delta)];
% TLable={'\fontname{宋体}墙体\fontname{Times new Roman}1';'2';'3';'4';'5';'6';'7';'8';'9';'\fontname{宋体}墙体\fontname{Times new Roman}10';'11';'12';'13';'14';'15';'16';'17';'18';'19'};
% set(gca,'FontName','Times new Roman','FontSize',16,'ThetaTick',TTick,'ThetaTickLabel',TLable,'Rlim',[0,1],'Rtick',[0:0.2:1.0],'RAxisLocation',45,'RTickLabel',{});
% L1=legend('','','','');
% set(L1,'FontSize',15);
% legend('boxoff');
% set(gcf,'units','normalized','position',[0.3 0.1 0.35 0.5]);
% print('-f1','H:\OneDrive\MathWork\Figures\f','-painters','-dmeta','-r600');

%% 动态演示1 GIF 
axis([-2,2,-2,2]);xlabel('X轴');ylabel('Y轴');box on;
t=0:0.02:10;  
Nt=size(t,2);
x=2*cos(t(1:Nt));
y=sin(t(1:Nt));
for i=1:Nt;
    plot(x(i),y(i),'o');
    hold on;
    %drawnow  %现在就画
    %pause(0.01)  %画完当前，停留0.01秒
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if i==1
         imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf,'DelayTime',1e-4);
    else
         imwrite(imind,cm,'test.gif','gif','WriteMode','append','DelayTime',1e-4);
    end
end

%% 动态演示2 movie
x=-8:0.5:8;
[XX,YY]=meshgrid(x);
r=sqrt(XX.^2+YY.^2)+eps;
Z=sin(r)./r;
surf(Z);
theAxes=axis;
fmat=moviein(20);
for j=1:20;
  surf(sin(2*pi*j/20)*Z,Z);
  axis(theAxes);
  fmat(:,j)=getframe;
end
movie(fmat,10)
