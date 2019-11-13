clear all;clc;close all;

%%%%%%%%%%%data prepare
startup;
Emeanf=0.96;
Edf=sqrt(0.06);
meanf=exp(Emeanf+0.5*Edf*Edf);
df=exp(2*Emeanf+Edf*Edf)*(exp(Edf*Edf)-1);
meanc=0.42*3.6;
dc=0.1;
DeltaT=0.26;%实验结果变异系数
A=326400;
Vk=10610;
NumSpl=300;
% Vmiu=1.06*Vk;
% Valpha=(pi/sqrt(6))/(Vmiu*0.3);
% Vbeta=Vmiu-0.5772/Valpha;


% N=100;
% f2s=rand(N,1);
% f2=logninv(f2s,Emeanf,Edf);
% cs=rand(N,1);
% c=norminv(cs,meanc/3.6,dc);
% X=[f2 c];
% R=0.1*sqrt(X(:,1)+sqrt(X(:,1)).*X(:,2));

X=[lognrnd(Emeanf,Edf,[NumSpl,1]) normrnd(meanc,dc,[NumSpl,1])];
Rtemp=0.1*sqrt(X(:,1)+sqrt(X(:,1)).*X(:,2));
R=Rtemp-mean(Rtemp)*DeltaT+mean(Rtemp)*DeltaT*2*rand(NumSpl,1);
xs=[meanf meanc];

% Xtrain=[2.31 	0.70 
% 4.67 	0.60 
% 2.07 	0.60 
% 3.69 	0.15 
% 8.76 	0.50 
% 2.82 	1.20 
% 4.67 	1.20 
% 1.08 	0.64 
% 0.52 	0.35 
% 3.39 	0.70 
% 3.39 	1.30 
% 0.10 	0.71 
% 0.77 	0.96 
% 1.08 	0.15 
% 1.08 	0.30 
% 1.08 	0.80 
% 2.07 	1.00 
% 5.02 	0.60 
% 5.02 	0.30 
% 0.31 	0.46 
% 0.92 	0.46 
% 5.76 	0.90 
% 5.76 	0.22 
% 0.23 	0.06 
% 0.23 	0.30 
% 0.77 	0.06 
% 0.77 	0.30 
% 0.77 	0.20 
% 0.64 	0.06 
% 0.64 	0.30 
% 0.64 	0.20 
% 1.08 	0.06 
% 1.08 	0.46 
% 1.08 	0.30 
% 1.08 	0.22 
% 1.08 	0.28 
% 1.08 	0.19 
% 1.08 	0.62 
% 1.08 	0.65 
% 1.08 	0.65 
% 1.08 	1.04 
% 1.08 	1.04 
% 1.08 	1.02
% 4.69    0.43
% 7.49    0.54
% 7.72    0.55
% 8.42    0.54
% 8.55    0.52
% 4.37    0.42
% ];
% Rtrain=[0.49 
% 0.66 
% 0.59 
% 0.44 
% 0.64 
% 0.31 
% 0.67 
% 0.29 
% 0.24 
% 0.48 
% 0.44 
% 0.32 
% 0.37 
% 0.14 
% 0.19 
% 0.65 
% 0.73 
% 0.52 
% 0.45 
% 0.32 
% 0.38 
% 0.76 
% 0.41 
% 0.06 
% 0.16 
% 0.07 
% 0.19 
% 0.16 
% 0.09 
% 0.21 
% 0.16 
% 0.11 
% 0.35 
% 0.28 
% 0.31 
% 0.42 
% 0.30 
% 0.66 
% 0.68 
% 0.73 
% 0.99 
% 0.95 
% 0.82 
% 0.74
% 0.97
% 0.75
% 0.67
% 0.81
% 0.38
% ];
% X=[X];
% R=[R];


%%%%%%%  training inputs
%prior models
        mL = {@meanLinear}; hypL = [1;1];

        M1=[true,false];
        M2=[false,true];
         
        mm1 = {'meanMask',M1,mL};   hypma1 = hypL(M1);
        mpo = {'meanPow',0.5,mm1};   hyppo = hypma1;

        mm2 = {'meanMask',M2,mL};   hypma2 = hypL(M2);
        mpr = {@meanProd,{mpo,mm2}}; hyppr = [hyppo; hypma2];
        
        msu = {'meanSum',{mm1,mpr}};  hypsu = [hypma1; hyppr];
        mpor = {'meanPow',0.5,msu};   hyppor = hypsu;

        
%       M3=[false,false,true,];
%       mm3={'meanMask',M3,mL};   hypma3 = hypL(M3);
%       msuF = {'meanSum',{mprod,mm3}};  hypsuF = [hypprod;  hypma3];
%       meanfunc = msuF;  
%       hyp.mean=hypsuF;
      
%         meanfunc = [];              % do not use a mean function
%         hyp.mean= [];
        meanfunc=mpor;
        hyp.mean= hyppor;
        covfunc = @covSEiso;            % Squared Exponental covariance function
        hyp.cov=[0 0];
        likfunc = @likGauss;            % Gaussian likelihood
        hyp.lik=0;    
         

%%%posterior hyperparameters
hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, X, R);
disp(hyp2); 

%%%%%%%%%%%%%Posterior Distributions of G
[Ymu Ys2 Fmu Fs2]= gp(hyp2, @infGaussLik, meanfunc, covfunc,likfunc,X,R,xs);
Rmu=A*Ymu;
Rs2=A*A*Ys2;

%%%通过多重积分求可靠度指标
% syms x y;
% RR=1-int(((2*pi*Rs2).^(-0.5))*exp(-1*(x-Rmu).^2/(2*Rs2)),x,-inf,y); 
% SS=int(RR*(2.35/(0.385*Vk))*(y/(0.385*Vk)).^(-3.35)).*exp(-1*(y/(0.385*Vk)).^(-2.35),y,0,inf);
% Pf=vpa(SS,2);

%%%通过卷积求可靠度指标
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
% %%%通过Fourier变换求可靠度指标
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

%%%%%%%%%%%%%plot 2-D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%export data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
print('-f1','H:\OneDrive\MathWork\Figures\f6','-painters','-dmeta','-r600');
% print('-f2','H:\MathWork\matlab\f2','-painters','-dmeta','-r600');






