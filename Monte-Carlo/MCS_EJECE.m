clear all;clc;close all;

%% input data
N=1e7;
Emeanf=1.43;
Edf=sqrt(0.12);
% meanf=exp(Emeanf+0.5*Edf*Edf);
% df=exp(2*Emeanf+Edf*Edf)*(exp(Edf*Edf)-1);
meanc=0.52;
dc=0.01;
A=477600;
Vk=18340;
DeltaT=0.26;%墙体名义抗剪强度实验结果变异系数

%% Monte Carlo Simuluation
f2=lognrnd(Emeanf,Edf,[N,1]);
c=normrnd(meanc,dc,[N,1]);
V=gevrnd(1/2.35,0.385*Vk/2.35,0.385*Vk,[N,1]);
Rtemp=0.1*A*sqrt(f2+3.6*c.*sqrt(f2));
R=Rtemp-mean(Rtemp)*DeltaT+mean(Rtemp)*DeltaT*2*rand(N,1);
H=(R<V);
Pf=sum(H)/N;
Beta=norminv(1-Pf);

%% 三维极限状态曲面图
% %当量正态化f2
% SS=makedist('Lognormal','mu',Emeanf,'sigma',Edf);
% SigmaN=normpdf(norminv(cdf(SS,Emeanf),0,1))/pdf(SS,Emeanf);
% miuN=Emeanf-SigmaN*norminv(cdf(SS,Emeanf),0,1);
% [Xcor Ycor]=meshgrid(linspace(3,4,100),linspace(0.6,0.7,100));
% Xcor=(Xcor-miuN)/SigmaN;
% Ycor=(Ycor-meanc)/dc;
% V=0.1*A*sqrt(Xcor+3.6*Ycor.*sqrt(Xcor));
% grid on;
% % for i=1:36
% %    camorbit(10,0);
% %    drawnow;
% % end

%% 三维柱状图
% %当量正态化V
%  SS=makedist('GeneralizedExtremeValue', 'k',1/2.35,'sigma',0.385*Vk/2.35,'mu',0.385*Vk);
%  SigmaV=normpdf(norminv(cdf(SS,mean(SS)),0,1))/pdf(SS,mean(SS));
%  miuV=Emeanf-SigmaV*norminv(cdf(SS,mean(SS)),0,1);
%  SV=normrnd(miuV,SigmaV,[N,1]);
% dat=[(R-mean(R))/std(R) (SV-miuV)/SigmaV];
% nBar=60;
% hist3(dat,[nBar nBar]);
% t=hist3(dat,[nBar nBar]);
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
% % % hold on;
% % % n = hist3(dat,[nBar nBar]); % default is to 10x10 bins
% % % n1 = n';
% % % n1(size(n,1) + 1, size(n,2) + 1) = 0;
% % % xb = linspace(mi
% % % n(dat(:,1)),max(dat(:,1)),size(n,1)+1);
% % % yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
% % % h = pcolor(xb,yb,n1);
% % % h.ZData = ones(size(n1)) * -max(max(n));
% % colormap(jet) % color map
% title('Data Point Density Histogram and Intensity Map');
% grid on;
% colorbar;
% view(3)

%% 平面柱状图
figure(1);

SigmaR=sqrt(217130000);
MiuR=95610;
RR=makedist('Normal',MiuR,SigmaR);
RoS=linspace(0,200000,1000);
h1=plot(RoS,pdf(RR,RoS),'LineWidth',1.5);
hold on;

nBar=100;
Xedges=[min(Rtemp):(max(Rtemp)-min(Rtemp))/nBar:max(Rtemp)]';
h4=histogram(Rtemp,Xedges,'Normalization','pdf','FaceColor','y');

set(gca,'FontName','Times','FontSize',13, 'xlim',[0,2.5e5],'xtick',[0:0.5e5:2.5e5],'ylim',[0,35e-6],'ytick',[0:5e-6:35e-6]);
box off;
%   xlabel('\itR\rm(kN)','FontName','Timesnewroman');
%   ylabel('\itS\rm(kN)','FontName','Timesnewroman');
%   zlabel('{\pi}_{\itR}、{\pi}_{\itS}','FontName','Timesnewroman');

  %修改坐标轴刻度显示方式 
  x_formatstring = '%5.0f';
  xtick = get(gca, 'xtick');
  y_formatstring = '%5.0d';
  ytick = get(gca, 'ytick');
  
  for i = 2:length(xtick)
    xticklabel{i} = sprintf(x_formatstring, xtick(i)/1000);
  end
  for j = 2:length(ytick)
    yticklabel{j} = sprintf(y_formatstring, ytick(j)*1e4);
  end
 set(gca, 'xticklabel', xticklabel); 
 set(gcf,'units','centimeters','position',[0 0 7 6]);

%% 散点图
% scatter3(f2,c,log(V/A));
% set(gca,'FontName','Timesnewroman','xlim',[0,max(f2)],'xtick',[0:max(f2)/5:max(f2)],'ylim',[0,max(c)],'ytick',[0:max(c)/5:max(c)],'zlim',[0,max(log(V/A))],'ztick',[0:max(log(V/A))/5:max(log(V/A))]);

%% export data
  print('-f1','C:\Users\Bin Peng\OneDrive\MathWork\Figures\f4(a)','-painters','-dmeta','-r600');