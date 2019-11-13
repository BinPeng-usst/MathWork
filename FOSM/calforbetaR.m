function [ alpha1 , xopt1 ,result,inum] = calforbetaR( initvalue , beta0 , alpha0 , xopt0 , jacg , g, iter,Cov,PdofS) %对选取的beta进行计算
result = 0; 
alpha = alpha0; 
xopt = xopt0; 
%initvalue为初始值 
miu = initvalue(1,:);%第一行为均值 
v = initvalue(2,:);%第二行为变异系数 
sgma = initvalue(3,:);%第三行为方差 
inum=iter; 
eps = 0.1;


while 1 
    inum=inum+1;
 %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%当量正态化%%%%%%%%%%%%%%%%%%%%%%%%%  
%  miu(1,1)=RF_Miu(xopt(1,1),'Lognormal1');%当量正态化的均值
%  miu(1,2)=RF_Miu(xopt(1,2),'Lognormal2');%当量正态化的均值
   miu(1,3)=RF_Miu(xopt(1,3),PdofS);%当量正态化的均值
%  sgma(1,1)=RF_Sigma(xopt(1,1),'Lognormal1');%当量正态化的标准差
%  sgma(1,2)=RF_Sigma(xopt(1,2),'Lognormal2');%当量正态化的标准差
   sgma(1,3)=RF_Sigma(xopt(1,3),PdofS);%当量正态化的标准差
%  Cov=Nataf_Cov(Cov,'Normal','Lognormal');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%当量正态化%%%%%%%%%%%%%%%%%%%%%%%%% 
%计算alpha和新的验算点xopt 
   Ggrad = jacgfunc(jacg,xopt); 
   alpha0 = (sgma .* Ggrad) / sqrt((sgma.*Ggrad)*Cov*(sgma.*Ggrad)'); 
   xopt0=miu-beta0 *sgma.*(Cov*alpha0')';
   sum1 = norm((alpha - alpha0),2); 
   sum2 = norm((xopt - xopt0),2); 

%验算点xi达到精度则退出循环 
   if ((sum1<eps*norm(alpha0,2))&&(sum2<eps*norm(xopt0,2)))
      alpha1 = alpha0; 
      xopt1 = xopt0; 
        if abs(gfunc(g,xopt0)) < eps 
          result = 1;  %功能函数同时达到精度,已计算出可靠指标 
         end 
      break; 
   end
   
   alpha = alpha0;
   xopt = xopt0;
end
