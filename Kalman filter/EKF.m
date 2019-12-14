clear all;clc;close all;
%%Preparation of the ambinet test data.IO is the detrended&truncated original signal, NS is SVDed signal, FS is EKFed signal
% Read data  
ImDataFile='C:\Users\pengbin\OneDrive - usst.edu.cn\桌面\Publication\振动与冲击\wall data\W1_data.mat';
DRepos=importdata(ImDataFile);
% Fst=1;Lnth=4095;Lst=Fst+Lnth;SclFtr=0.01;
% IOA=[DRepos.A_input1T(Fst:Lst),DRepos.A_output1T(Fst:Lst),DRepos.A_output2T(Fst:Lst),DRepos.A_output3T(Fst:Lst),DRepos.A_output4T(Fst:Lst),DRepos.A_output5T(Fst:Lst),DRepos.A_output6T(Fst:Lst),DRepos.A_output7T(Fst:Lst)];
% IOB=[DRepos.B_input1T(Fst:Lst),DRepos.B_output1T(Fst:Lst),DRepos.B_output2T(Fst:Lst),DRepos.B_output3T(Fst:Lst),DRepos.B_output4T(Fst:Lst),DRepos.B_output5T(Fst:Lst),DRepos.B_output6T(Fst:Lst),DRepos.B_output7T(Fst:Lst)];
% IO1=[DRepos.C_input1T(Fst:Lst),DRepos.C_output1T(Fst:Lst),DRepos.C_output2T(Fst:Lst),DRepos.C_output3T(Fst:Lst),DRepos.C_output4T(Fst:Lst),DRepos.C_output5T(Fst:Lst),DRepos.C_output6T(Fst:Lst),DRepos.C_output7T(Fst:Lst)];
% IO2=[DRepos.D_input1T(Fst:Lst),DRepos.D_output1T(Fst:Lst),DRepos.D_output2T(Fst:Lst),DRepos.D_output3T(Fst:Lst),DRepos.D_output4T(Fst:Lst),DRepos.D_output5T(Fst:Lst),DRepos.D_output6T(Fst:Lst),DRepos.D_output7T(Fst:Lst)];
% IO3=[DRepos.E_input1T(Fst:Lst),DRepos.E_output1T(Fst:Lst),DRepos.E_output2T(Fst:Lst),DRepos.E_output3T(Fst:Lst),DRepos.E_output4T(Fst:Lst),DRepos.E_output5T(Fst:Lst),DRepos.E_output6T(Fst:Lst),DRepos.E_output7T(Fst:Lst)];
Fst=1;Lnth=size(DRepos.A_input,1)-1;Lst=Fst+Lnth;SclFtr=0.01;
IOA=[DRepos.A_input(Fst:Lst),DRepos.A_output1(Fst:Lst),DRepos.A_output2(Fst:Lst),DRepos.A_output3(Fst:Lst),DRepos.A_output4(Fst:Lst),DRepos.A_output5(Fst:Lst),DRepos.A_output6(Fst:Lst),DRepos.A_output7(Fst:Lst)];

Fst=1;Lnth=size(DRepos.B_input,1)-1;Lst=Fst+Lnth;SclFtr=0.01;
IOB=[DRepos.B_input(Fst:Lst),DRepos.B_output1(Fst:Lst),DRepos.B_output2(Fst:Lst),DRepos.B_output3(Fst:Lst),DRepos.B_output4(Fst:Lst),DRepos.B_output5(Fst:Lst),DRepos.B_output6(Fst:Lst),DRepos.B_output7(Fst:Lst)];

Fst=1;Lnth=size(DRepos.C_input,1)-1;Lst=Fst+Lnth;SclFtr=0.01;
IO1=[DRepos.C_input(Fst:Lst),DRepos.C_output1(Fst:Lst),DRepos.C_output2(Fst:Lst),DRepos.C_output3(Fst:Lst),DRepos.C_output4(Fst:Lst),DRepos.C_output5(Fst:Lst),DRepos.C_output6(Fst:Lst),DRepos.C_output7(Fst:Lst)];

Fst=1;Lnth=size(DRepos.D_input,1)-1;Lst=Fst+Lnth;SclFtr=0.01;
IO2=[DRepos.D_input(Fst:Lst),DRepos.D_output1(Fst:Lst),DRepos.D_output2(Fst:Lst),DRepos.D_output3(Fst:Lst),DRepos.D_output4(Fst:Lst),DRepos.D_output5(Fst:Lst),DRepos.D_output6(Fst:Lst),DRepos.D_output7(Fst:Lst)];

Fst=1;Lnth=size(DRepos.E_input,1)-1;Lst=Fst+Lnth;SclFtr=0.01;
IO3=[DRepos.E_input(Fst:Lst),DRepos.E_output1(Fst:Lst),DRepos.E_output2(Fst:Lst),DRepos.E_output3(Fst:Lst),DRepos.E_output4(Fst:Lst),DRepos.E_output5(Fst:Lst),DRepos.E_output6(Fst:Lst),DRepos.E_output7(Fst:Lst)];

IO=IOA*SclFtr;
MsmtNse=mean(std(IO))/size(IO,2);
PcsNse=0.1*MsmtNse;

% Specify ambient test parameters
SplFreqcy=200;
DeltaT=1/SplFreqcy;
t=[0:1/SplFreqcy:(size(IO,1)-1)/SplFreqcy]';
% SenseOfAcc=32.8 ;%Gal/v
% SenseOfInstru=34.2; %Gal/v
% GainOfInstru=1;
% R=1/((SenseOfInstru/SenseOfAcc)*GainOfInstru);

% Detrend
DtrdOdr=2;
for i=1:size(IO,2)
    IO(:,i)=dtrend(IO(:,i)-polyval(polyfit(t,IO(:,i),DtrdOdr),t));
end

% Delete small voltage turbulents
% for i=1:size(IO,2)
%     for j=1:size(IO(:,i))
%       if abs(IO(j,i))>1e-3
%         IO(j,i)=IO(j,i);
%         else
%          IO(j,i)=0;
%         end
%     end
% end

% Integrating the acceleration to get velocity and dislacement
for i=1:size(IO,2)
V_th(:,i)=cumtrapz(t,IO(:,i));
end

for j=1:size(V_th,2)
D_th(:,j)=cumtrapz(t,V_th(:,j));
end

% Form the measurement matrix
[U1,e1,V1] = svd(D_th(:,2:8),0);
NS1=U1(:,1)*e1(1,1)*V1(1,:);
[U2,e2,V2] = svd(V_th(:,2:8),0);
NS2=U2(:,1)*e2(1,1)*V2(1,:);
NS=[NS1,NS2];
% % SVD decompositon the measurement matrix
% [U,e,V] = svd(NS,0);
% NS=U(:,1)*e(1,1)*V(1,:);

%% Preparation of structural matricies
M=DRepos.M; 
C=DRepos.C;
K=DRepos.K; 
D=[DRepos.juzhen;DRepos.juzhen]; 

%% EKF construction

Phi=[zeros(size(M)),eye(size(M));(-1)*(M^-1)*K,(-1)*(M^-1)*C];
Psi=[zeros(size(M));(-1)*M];

A=expm(Phi*DeltaT);
B=(Phi^-1)*(1-(expm(Phi*DeltaT))^-1)*Psi;

StateFcn=@(X,U)(A*X+B*U);
MeasurementFcn=@(Y)(D*Y);
obj = extendedKalmanFilter(StateFcn,MeasurementFcn,zeros(size(A,2),1),'StateCovariance',5);
obj.ProcessNoise=PcsNse;
obj.MeasurementNoise=MsmtNse; 

for k = 1:size(NS,1)
  [CorrectedState,CorrectedStateCovariance] = correct(obj,NS(k,1:14)'); 
  [PredictedState,PredictedStateCovariance] = predict(obj,IO(k,1)*ones(size(B,2),1));
  FS(k,1:14)=(D*CorrectedState)';
end

%% PSD of the filtered signal
[PsdT_1,f0]=pwelch(D_th(:,2),[],[],[],SplFreqcy);%Hamming窗，默认窗长度、重叠长度和DFT点数
[PsdN_1,f1]=pwelch(NS(:,1),[],[],[],SplFreqcy);%Hamming窗，默认窗长度、重叠长度和DFT点数
[PsdF_1,f2]=pwelch(FS(:,1),[],[],[],SplFreqcy);%Hamming窗，默认窗长度、重叠长度和DFT点数
% PsdT_1(1)=0;
% PsdF_1(1)=0;

%% Output
figure(1)
subplot(2,2,[1,2]);
for i=1:8
    plot3(t,i*ones(size(t,1),1),D_th(:,i));hold on;
end
    
subplot(2,2,3);
plot(t,D_th(:,2)); hold on;
plot(t,NS(:,1)+max(NS(:,1)+max(D_th(:,1)))); 
plot(t,FS(:,1)+max(FS(:,1)+max(NS(:,1))+max(D_th(:,1)))); 
legend('\fontname{宋体}预处理后原纪录\fontname{Times new Roman}(D-th)','SVD\fontname{宋体}分解重构\fontname{Times new Roman}(NS)','EKF\fontname{宋体}滤波后\fontname{Times new Roman}(FS)','Location','best');
xlabel('\fontname{宋体}时间\fontname{Times new Roman}(s)','FontSize',10);
ylabel('\fontname{宋体}变形\fontname{Times new Roman}(m)','FontSize',10);
% set(gca,'FontName','Times new Roman','FontSize',11);m
% set(gcf,'Units','centimeters','Position',[0 0 16 16],'Resize','off');
% RFile='C:\Users\pengbin\Desktop\T';
% print('-f1',RFile,'-painters','-dmeta','-r600');
subplot(2,2,4);
[pks0,loc0]=max(PsdT_1);BscFreq0=f0(loc0);
[pks1,loc1]=max(PsdN_1);BscFreq1=f1(loc1);
[pks2,loc2]=max(PsdF_1);BscFreq2=f2(loc2);
plot(f0,PsdT_1);hold on;
plot(f1,PsdN_1);hold on;
plot(f2,PsdF_1);hold on;
legend(['PSD of D-th (Hz)',' ',num2str(BscFreq0)],['PSD of NS (Hz)',' ',num2str(BscFreq1)],['PSD of FS (Hz)',' ',num2str(BscFreq2)],'Location','best');