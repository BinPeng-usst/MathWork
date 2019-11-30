clear all;clc;close all;
%%Preparation of the ambinet test data.S is the original signal,TIO is the detrended&truncated signal,NS is SVD-processed signal,FS is EKF filtered signal
% Read data  
ImDataFile='C:\Users\pengbin\OneDrive - usst.edu.cn\桌面\Publication\振动与冲击\wall data\W1_data.mat';
DRepos=importdata(ImDataFile);
IO=[DRepos.A_input(115200:115200+65536),DRepos.A_output1(115200:115200+65536),DRepos.A_output2(115200:115200+65536),DRepos.A_output3(115200:115200+65536),DRepos.A_output4(115200:115200+65536),DRepos.A_output5(115200:115200+65536),DRepos.A_output6(115200:115200+65536),DRepos.A_output7(115200:115200+65536)];

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
for i=1:size(IO,2)
    for j=1:size(IO(:,i))
      if abs(IO(j,i))>1e-3
        TIO(j,i)=IO(j,i);
        else
         TIO(j,i)=0;
        end
    end
end

% Integrating the acceleration to get velocity and dislacement
for i=1:size(TIO,2)
 V_th(:,i)=cumtrapz(t,TIO(:,i));
end

for j=1:size(V_th,2)
 D_th(:,j)=cumtrapz(t,V_th(:,j));
end

% SVD decompositon the measured records
NS=[D_th(:,2:8),V_th(:,2:8)];
[U,e,V] = svd(NS,0);
NS=U(:,1)*e(1,1)*V(1,:);

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
obj = extendedKalmanFilter(StateFcn,MeasurementFcn,NS(1,1)*ones(size(A,2),1),'StateCovariance',5);
obj.ProcessNoise=0.618;
obj.MeasurementNoise=1; 

for k = 1:size(NS,1)
  [CorrectedState,CorrectedStateCovariance] = correct(obj,NS(k,1:14)'); 
  [PredictedState,PredictedStateCovariance] = predict(obj,TIO(k,1)*ones(size(B,2),1));
  FS(k,1:14)=(D*CorrectedState)';
end

%% PSD of the filtered signal
[PsdT_1,f0]=pwelch(TIO(:,2),[],[],1024,SplFreqcy);%Hamming窗，默认窗长度、重叠长度和DFT点数
[PsdF_1,f1]=pwelch(FS(:,1),[],[],1024,SplFreqcy);%Hamming窗，默认窗长度、重叠长度和DFT点数
PsdT_1(1:20)=0;
PsdF_1(1:20)=0;

%% Output
figure(1)
plot(t,TIO(:,2)); hold on;
plot(t,FS(:,1)); 
legend('预处理后原纪录（TIO）','EKF滤波后（FS）');
xlabel('\fontname{宋体}时间\fontname{Times new Roman}(s)','FontSize',10);
ylabel('\fontname{宋体}加速度\fontname{Times new Roman}(cm/s^{2})','FontSize',10);

% ax=gca;outerpos=ax.OuterPosition;ti=ax.TightInset; 
% left=outerpos(1)+ti(1);bottom=outerpos(2)+ti(2);ax_width=outerpos(3)-ti(1)-ti(3);ax_height=outerpos(4)-ti(2)-ti(4);
% ax.Position=[left bottom ax_width ax_height];
set(gca,'FontName','Times new Roman','FontSize',11);
set(gcf,'Units','centimeters','Position',[0 0 16 16],'Resize','off');
RFile='C:\Users\Bin Peng\Desktop\T';
print('-f1',RFile,'-painters','-dmeta','-r600');

figure(2);
[pks0,loc0]=max(PsdT_1);BsFreq0=f0(loc0);
[pks1,loc1]=max(PsdF_1);BsFreq1=f1(loc1);
plot(f0,PsdT_1);hold on;
plot(f1,PsdF_1);hold on;
legend(['PSD of TIO',' ',num2str(BsFreq0)],['PSD of FS',' ',num2str(BsFreq2)]);


