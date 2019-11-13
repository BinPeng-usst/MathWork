%����ΪstrRlbt.m,��Matlab�����������strRlbtִ�м��� 
clear all;clc;close all;
miu = [2.70 0.42 0.597*10610];%��ֵ 
v = [0.26 0.24 0.3];%����ϵ�� 
sgma = miu .* v;%���� 
Cov=[1.0 0.0 0.0;%���ϵ�� 
     0.0 1.0 0.0;
     0.0 0.0 1.0];
eps = 0.1; %���� 

Vk=10610;
PdofS=makedist('GeneralizedExtremeValue', 'k',1/2.35,'sigma',0.385*Vk/2.35,'mu',0.385*Vk);%����ЧӦ�ķֲ�
%PdofS=makedist('Lognormal', 'mu',3.681,'sigma',0.125);
%PdofS=makedist('Lognormal', 'mu',3.910,'sigma',0.05);
% PdofS=makedist('ExtremeValue', 'mu',-910.0,'sigma',1/0.006413);
      
syms f c V; 
g = sym('0.1*326400*sqrt(f+3.6*c*sqrt(f))-V');%���ܺ��� 
jacg = jacobian(g ,[f;c;V]);%�����ſɱȾ��� 
initvalue = [miu;v;sgma];%������������ 
Tnum=0; 


%ѡȡbeta������x0=miu 
beta1 = 3;
xopt = miu; 
alpha = zeros(1,3); 
[alpha1 , xopt1 , result, inum] = calforbetaR(initvalue , beta1 , alpha, xopt , jacg , g, 0, Cov,PdofS); % ����������ʱ
Tnum=Tnum+inum;
if result == 1 
   disp('�ɿ�ָ��Ϊ��һ�μٶ���ֵ'); 
   disp(beta1);
   disp('���������Ϊ��'); 
   disp(xopt1); 
   disp('������㴦���ܺ���ֵΪ��'); 
   disp(gfunc(g,xopt1));  
   disp('�ܵ�������Ϊ'); 
   disp(Tnum);
   return 
end 
 
%�ٴμٶ�beta 
beta2 = 2.5; 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%������̬��%%%%%%%%%%%%%%%%%%%%%%%%%
%  miu(1,1)=RF_Miu(xopt1(1,1),'Lognormal1');%������̬���ľ�ֵ
%  miu(1,2)=RF_Miu(xopt1(1,2),'Lognormal2');%������̬���ľ�ֵ
   miu(1,3)=RF_Miu(xopt1(1,3),PdofS);%������̬���ľ�ֵ
%  sgma(1,1)=RF_Sigma(xopt1(1,1),'Lognormal1');%������̬���ı�׼��
%  sgma(1,2)=RF_Sigma(xopt1(1,2),'Lognormal2');%������̬���ı�׼��
   sgma(1,3)=RF_Sigma(xopt1(1,3),PdofS);%������̬���ı�׼��
%  Cov=Nataf_Cov(Cov,'Normal','Lognormal');
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%������̬��%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Ggrad = jacgfunc(jacg,xopt1);
   alpha = (sgma .* Ggrad) / sqrt((sgma.*Ggrad)*Cov*(sgma.*Ggrad)'); 
   xopt=miu-beta2 *sgma.*(Cov*alpha')';
  [alpha2 , xopt2 , result, inum] = calforbetaR( initvalue , beta2 , alpha , xopt , jacg , g, inum, Cov,PdofS); 
  Tnum=Tnum+inum;
if result == 1 
   disp('�ɿ�ָ��Ϊ�ڶ��μٶ���ֵ'); 
   disp(beta2);
   disp('���������Ϊ��'); 
   disp(xopt2); 
   disp('������㴦���ܺ���ֵΪ��'); 
   disp(gfunc(g,xopt2)); 
   disp('��������Ϊ'); 
   disp(Tnum); 
   return
end

%beta�������
g1 = gfunc(g,xopt1); 
g2 = gfunc(g,xopt2); 
while ((abs(g2)> eps) | (beta2<0))% ���Կɿ���ָ�����
    temp = beta2; 
    beta2 = beta2 - (beta2 - beta1)/(g2 - g1) * g2; %���߷�
    beta1 = temp; 
   [alpha2 , xopt2 , result, inum] = calforbetaR( initvalue , beta2 , alpha2 , xopt2 , jacg , g, inum, Cov,PdofS);%�ڲ�������������ȵ���
     g1 = g2; 
   g2 = gfunc(g,xopt2); 
   Tnum=Tnum+inum;
   if (result == 1)&&(beta2>=0)
      break 
   end 
end 
 
disp('�ɿ�ָ��Ϊ��'); 
disp(beta2); 
disp('���������Ϊ��'); 
disp(xopt2); 
disp('������㴦���ܺ���ֵΪ��'); 
disp(g2);  
disp('�ܵ�������Ϊ'); 
disp(Tnum);  

