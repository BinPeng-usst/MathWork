function [ alpha1 , xopt1 ,result,inum] = calforbetaR( initvalue , beta0 , alpha0 , xopt0 , jacg , g, iter,Cov,PdofS) %��ѡȡ��beta���м���
result = 0; 
alpha = alpha0; 
xopt = xopt0; 
%initvalueΪ��ʼֵ 
miu = initvalue(1,:);%��һ��Ϊ��ֵ 
v = initvalue(2,:);%�ڶ���Ϊ����ϵ�� 
sgma = initvalue(3,:);%������Ϊ���� 
inum=iter; 
eps = 0.1;


while 1 
    inum=inum+1;
 %  %%%%%%%%%%%%%%%%%%%%%%%%%%%%������̬��%%%%%%%%%%%%%%%%%%%%%%%%%  
%  miu(1,1)=RF_Miu(xopt(1,1),'Lognormal1');%������̬���ľ�ֵ
%  miu(1,2)=RF_Miu(xopt(1,2),'Lognormal2');%������̬���ľ�ֵ
   miu(1,3)=RF_Miu(xopt(1,3),PdofS);%������̬���ľ�ֵ
%  sgma(1,1)=RF_Sigma(xopt(1,1),'Lognormal1');%������̬���ı�׼��
%  sgma(1,2)=RF_Sigma(xopt(1,2),'Lognormal2');%������̬���ı�׼��
   sgma(1,3)=RF_Sigma(xopt(1,3),PdofS);%������̬���ı�׼��
%  Cov=Nataf_Cov(Cov,'Normal','Lognormal');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%������̬��%%%%%%%%%%%%%%%%%%%%%%%%% 
%����alpha���µ������xopt 
   Ggrad = jacgfunc(jacg,xopt); 
   alpha0 = (sgma .* Ggrad) / sqrt((sgma.*Ggrad)*Cov*(sgma.*Ggrad)'); 
   xopt0=miu-beta0 *sgma.*(Cov*alpha0')';
   sum1 = norm((alpha - alpha0),2); 
   sum2 = norm((xopt - xopt0),2); 

%�����xi�ﵽ�������˳�ѭ�� 
   if ((sum1<eps*norm(alpha0,2))&&(sum2<eps*norm(xopt0,2)))
      alpha1 = alpha0; 
      xopt1 = xopt0; 
        if abs(gfunc(g,xopt0)) < eps 
          result = 1;  %���ܺ���ͬʱ�ﵽ����,�Ѽ�����ɿ�ָ�� 
         end 
      break; 
   end
   
   alpha = alpha0;
   xopt = xopt0;
end
