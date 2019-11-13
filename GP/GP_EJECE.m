clear all;clc;close all;

%%
%data prepare
Emeanf=0.96;
Edf=sqrt(0.06);    
meanf=exp(Emeanf+0.5*Edf*Edf);
df=exp(2*Emeanf+Edf*Edf)*(exp(Edf*Edf)-1);
meanc=0.42;
dc=0.1;
DeltaT=0.26;%实验结果变异系数

A=326400;
Vk=10610;
NumSpl=300;
X=[lognrnd(Emeanf,Edf,[NumSpl,1]) normrnd(meanc,dc,[NumSpl,1])];
Rtemp=0.1*sqrt(X(:,1)+3.6*sqrt(X(:,1)).*X(:,2));
R=Rtemp-mean(Rtemp)*DeltaT+mean(Rtemp)*DeltaT*2*rand(NumSpl,1);
xs=[2.70 0.42];

%%
%%%%%  training inputs
%prior models
      Prior=@(Me) 0.1*sqrt(Me(:,1)+3.6*sqrt(Me(:,1)).*Me(:,2));           
      gprMdl = fitrgp(X,R,'BasisFunction',Prior,'Beta',[1],'KernelFunction','squaredexponential',...
      'Verbose',0);
      yres=resubPredict(gprMdl);
      [ypred,ysd,yint] = predict(gprMdl,xs,'Alpha',0.1);

%%%%%%%%%%%%%Posterior Distributions of G
Rmu=A*ypred;
Rs2=A^2*ysd^2;

%%
%%求可靠度指标
%%%多重积分法
% syms x y;
% RR=1-int(((2*pi*Rs2).^(-0.5))*exp(-1*(x-Rmu).^2/(2*Rs2)),x,-inf,y); 
% SS=int(RR*(2.35/(0.385*Vk))*(y/(0.385*Vk)).^(-3.35)).*exp(-1*(y/(0.385*Vk)).^(-2.35),y,0,inf);
% Pf=vpa(SS,2);

%%%卷积法
n=100000;
RoS=linspace(0,Rmu+10*sqrt(Rs2),n);
Delta=(Rmu+10*sqrt(Rs2))/n;
Pf=0;
RR=makedist('Normal',Rmu,sqrt(Rs2));
SS=makedist('GeneralizedExtremeValue', 'k',1/2.35,'sigma',0.385*Vk/2.35,'mu',0.385*Vk);
Rsample=pdf(RR,RoS);
Ssample=cdf(SS,RoS);
for i=1:n
    Pf=Pf+Rsample(i)*Delta*(1-Ssample(i));
end
Beta=norminv(1-Pf);

%%Fourier变换法
% syms x y;
% F1=fourier(((2*pi*Rs2).^(-0.5))*exp(-1*(x-Rmu).^2/(2*Rs2))); 
% F2=fourier((2.35/(0.385*Vk))*(y/(0.385*Vk)).^(-3.35)).*exp(-1*(y/(0.385*Vk)).^(-2.35));
% Pf=vpa(ifourier(F1*F2),2)
%%% 模拟指定分布的样本
% xS=linspace(0,Vk*10,100);
% pdf =
% proppdf = @(x,y) exp(-1*(0.5*x/(0.3*Vmiu))^(-2));
% proprnd = @(x)  unifrnd(xS(1),xS(10000));
% nsamples = 10000;
% [spl,accept]= mhsample(xS(1),nsamples,'pdf',pdf,'proprnd',proprnd,'proppdf',proppdf,'thin',21,'burnin',2000);
% [SS,xS] = ksdensity(spl,xS)
%%%求指定分布的交点
%  syms u
%  eqn = (1/sqrt(Rs2*2*pi))*exp(-(u-Rmu)*(u-Rmu)/(2*Rs2))==(2.35/(0.385*Vk))*((u/(0.385*Vk))^(-3.35))*exp(-1*(u/(0.385*Vk))^(-2.35));
%  lp=solve(eqn,u,'real', 1);
%%
%%%plot 2-D%%%
 figure(1);
 plot(RoS,pdf(RR,RoS),'LineWidth',5);
 hold on;
 plot(RoS,pdf(SS,RoS),'r','LineWidth',3,'linestyle',':');
 set(gca,'FontName','Timesnewroman','FontSize',35, 'xlim',[0,2.5e5],'xtick',[0:0.5e5:2.5e5],'ylim',[0,25e-5],'ytick',[0:5e-5:25e-5]);
 box off;
 %修改X轴刻度显示方式 
  x_formatstring = '%5.2d';
  xtick = get(gca, 'xtick');
  for i = 2:length(xtick)
    xticklabel{i} = sprintf(x_formatstring, xtick(i)/1000);
  end
  set(gca, 'xticklabel', xticklabel);    
  %修改X轴刻度显示方式
% legend({'Posterior PDF of the shear resistance R','PDF of the earthquake shear force Q'},'FontName','Timesnewroman','FontSize',21,'units','normalized','position',[0.46 0.7 0.2 0.2]);
 legend({'',''},'FontName','Timesnewroman','FontSize',21,'units','normalized','position',[0.26 0.6 0.2 0.2]);
 legend('boxoff');
%  ylabel('{\pi}_{\itR}、{\pi}_{\itS}','FontName','Timesnewroman');
%  xlabel('\itR、\itS\rm(kN)','FontName','Timesnewroman');
 set(gcf,'units','normalized','position',[0.3 0.1 0.39 0.5]);

% figure(2);
% f = [ymu+3*sqrt(ys2); flipdim(ymu-3*sqrt(ys2),1)];
% subplot(1,2,1);
% fill([X(:,1); flipdim(X(:,1),1)], f,'g','FaceAlpha',0.2,'EdgeColor','g','EdgeAlpha',0.2,'linestyle', '- .'); %
% subplot(1,2,2);
% fill([X(:,2); flipdim(X(:,2),1)], f,'g','FaceAlpha',0.2,'EdgeColor','g','EdgeAlpha',0.2,'linestyle', '- .');

%%
%%%%%%%%%%%%export data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('-f1','H:\OneDrive\MathWork\Figures\f6','-painters','-dmeta','-r600');
% print('-f2','H:\MathWork\matlab\f2','-painters','-dmeta','-r600');






