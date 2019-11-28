clear all;clc;close all;
%%Preparation of the ambinet test data.S is the original signal,TIO is the detrended&truncated signal,NS is SVD-processed signal,FS is EKF filtered signal
SamFreq=200;Rsl=1e-3;% Ambient test parameters
ImDataFile1='C:\Users\pengbin\OneDrive - usst.edu.cn\桌面\Publication\振动与冲击\data1.mat';
DRepos=importdata(ImDataFile1);
IO=[DRepos.input,DRepos.output1,DRepos.output2,DRepos.output3,DRepos.output4,DRepos.output5,DRepos.output6,DRepos.output7];

%detrend
DtrdOdr=2;
t=[0:1/SamFreq:(size(IO,1)-1)/SamFreq]';
for i=1:size(IO,2)
    IO(:,i)=dtrend(IO(:,i)-polyval(polyfit(t,IO(:,i),DtrdOdr),t));
end

%delete small voltage turbulents
for i=1:size(IO,2)
    for j=1:size(IO(:,i))
      if abs(IO(j,i))>Rsl
        TIO(j,i)=IO(j,i);
        j=j+1;
       else
         TIO(j,i)=0;
         j=j+1;
        end
    end
end

% SVD decompositon the output records
for j=2:8;
  NS(:,j-1)=TIO(:,j);
end
[U,e,V] = svd(NS,0);
NS=U(:,1)*e(1,1)*V(1,:);

%% Preparation of structural matricies
M=DRepos.M; 
C=DRepos.C;
K=DRepos.K; 
D=DRepos.juzhen; 

%% EKF construction
DeltaT=1/SamFreq;

Phi=[zeros(size(M)),eye(size(M));(-1)*(M^-1)*K,(-1)*(M^-1)*C];
Psi=[zeros(size(M));(-1)*M];

A=expm(Phi*DeltaT);
B=(Phi^-1)*(1-(expm(Phi*DeltaT))^-1)*Psi;

StateFcn=@(X,U)(A*X+B*U);
MeasurementFcn=@(Y)(D*Y);
obj = extendedKalmanFilter(StateFcn,MeasurementFcn,NS(1,1)*ones(size(A,2),1),'StateCovariance',5);
obj.ProcessNoise=0.618;
obj.MeasurementNoise=1; 

for k = 1:size(NS)
  [CorrectedState,CorrectedStateCovariance] = correct(obj,NS(k,1:7)'); 
  [PredictedState,PredictedStateCovariance] = predict(obj,TIO(k,1)*ones(size(B,2),1));
  FS(k,1:7)=(D*CorrectedState)';
end

%% PSD of the filtered signal
[PsdT_1,f0]=pwelch(TIO(:,2),[],[],1024,SamFreq);%Hamming窗，默认窗长度、重叠长度和DFT点数
[PsdN_1,f1]=pwelch(NS(:,1),[],[],1024,SamFreq);%Hamming窗，默认窗长度、重叠长度和DFT点数
[PsdF_1,f2]=pwelch(FS(:,1),[],[],1024,SamFreq);%Hamming窗，默认窗长度、重叠长度和DFT点数
PsdT_1(1:20)=0;
PsdN_1(1:20)=0;
PsdF_1(1:20)=0;

%% Output
figure(1)
plot(t,TIO(:,2)); hold on;
plot(t,NS(:,1)); hold on;
plot(t,FS(:,1)); 
legend('预处理后原纪录（TIO）','SVD分解后（NS）','EKF滤波后（FS）');
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
[pks1,loc1]=max(PsdN_1);BsFreq1=f1(loc1);
[pks2,loc2]=max(PsdF_1);BsFreq2=f2(loc2);
plot(f0,PsdT_1);hold on;
plot(f1,PsdN_1);hold on;
plot(f2,PsdF_1); 
legend(['PSD of TIO',' ',num2str(BsFreq0)],['PSD of NS',' ',num2str(BsFreq1)],['PSD of FS',' ',num2str(BsFreq2)]);
