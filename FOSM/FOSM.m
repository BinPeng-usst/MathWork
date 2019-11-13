%保存为strRlbt.m,在Matlab命令窗口中输入strRlbt执行即可 
clear all;clc;close all;
miu = [2.70 0.42 0.597*10610];%均值 
v = [0.26 0.24 0.3];%变异系数 
sgma = miu .* v;%方差 
Cov=[1.0 0.0 0.0;%相关系数 
     0.0 1.0 0.0;
     0.0 0.0 1.0];
eps = 0.1; %精度 

Vk=10610;
PdofS=makedist('GeneralizedExtremeValue', 'k',1/2.35,'sigma',0.385*Vk/2.35,'mu',0.385*Vk);%荷载效应的分布
%PdofS=makedist('Lognormal', 'mu',3.681,'sigma',0.125);
%PdofS=makedist('Lognormal', 'mu',3.910,'sigma',0.05);
% PdofS=makedist('ExtremeValue', 'mu',-910.0,'sigma',1/0.006413);
      
syms f c V; 
g = sym('0.1*326400*sqrt(f+3.6*c*sqrt(f))-V');%功能函数 
jacg = jacobian(g ,[f;c;V]);%计算雅可比矩阵 
initvalue = [miu;v;sgma];%用作函数参数 
Tnum=0; 


%选取beta定义x0=miu 
beta1 = 3;
xopt = miu; 
alpha = zeros(1,3); 
[alpha1 , xopt1 , result, inum] = calforbetaR(initvalue , beta1 , alpha, xopt , jacg , g, 0, Cov,PdofS); % 各量不独立时
Tnum=Tnum+inum;
if result == 1 
   disp('可靠指标为第一次假定的值'); 
   disp(beta1);
   disp('最终验算点为'); 
   disp(xopt1); 
   disp('在验算点处功能函数值为'); 
   disp(gfunc(g,xopt1));  
   disp('总迭代次数为'); 
   disp(Tnum);
   return 
end 
 
%再次假定beta 
beta2 = 2.5; 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%当量正态化%%%%%%%%%%%%%%%%%%%%%%%%%
%  miu(1,1)=RF_Miu(xopt1(1,1),'Lognormal1');%当量正态化的均值
%  miu(1,2)=RF_Miu(xopt1(1,2),'Lognormal2');%当量正态化的均值
   miu(1,3)=RF_Miu(xopt1(1,3),PdofS);%当量正态化的均值
%  sgma(1,1)=RF_Sigma(xopt1(1,1),'Lognormal1');%当量正态化的标准差
%  sgma(1,2)=RF_Sigma(xopt1(1,2),'Lognormal2');%当量正态化的标准差
   sgma(1,3)=RF_Sigma(xopt1(1,3),PdofS);%当量正态化的标准差
%  Cov=Nataf_Cov(Cov,'Normal','Lognormal');
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%当量正态化%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Ggrad = jacgfunc(jacg,xopt1);
   alpha = (sgma .* Ggrad) / sqrt((sgma.*Ggrad)*Cov*(sgma.*Ggrad)'); 
   xopt=miu-beta2 *sgma.*(Cov*alpha')';
  [alpha2 , xopt2 , result, inum] = calforbetaR( initvalue , beta2 , alpha , xopt , jacg , g, inum, Cov,PdofS); 
  Tnum=Tnum+inum;
if result == 1 
   disp('可靠指标为第二次假定的值'); 
   disp(beta2);
   disp('最终验算点为'); 
   disp(xopt2); 
   disp('在验算点处功能函数值为'); 
   disp(gfunc(g,xopt2)); 
   disp('迭代次数为'); 
   disp(Tnum); 
   return
end

%beta迭代求解
g1 = gfunc(g,xopt1); 
g2 = gfunc(g,xopt2); 
while ((abs(g2)> eps) | (beta2<0))% 外层对可靠度指标迭代
    temp = beta2; 
    beta2 = beta2 - (beta2 - beta1)/(g2 - g1) * g2; %切线法
    beta1 = temp; 
   [alpha2 , xopt2 , result, inum] = calforbetaR( initvalue , beta2 , alpha2 , xopt2 , jacg , g, inum, Cov,PdofS);%内层对验算点和灵敏度迭代
     g1 = g2; 
   g2 = gfunc(g,xopt2); 
   Tnum=Tnum+inum;
   if (result == 1)&&(beta2>=0)
      break 
   end 
end 
 
disp('可靠指标为'); 
disp(beta2); 
disp('最终验算点为'); 
disp(xopt2); 
disp('在验算点处功能函数值为'); 
disp(g2);  
disp('总迭代次数为'); 
disp(Tnum);  

