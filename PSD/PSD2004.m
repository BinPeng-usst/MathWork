clear all;clc;close all;
% %% Data input
% % Sampling Parameters
% % ImportSourceFile1='C:\Users\Administrator\Desktop\zhong\1��ͨ�� 3.dat';
% % ImportSourceFile2='C:\Users\Administrator\Desktop\zhong\1��ͨ�� 3.dat';
% % ImportSourceFile3='C:\Users\Administrator\Desktop\zhong\1��ͨ�� 3.dat';
% ImportResponseFile1='C:\Users\2013-FWQ\Desktop\AI1-01_20190513084355.txt';
% ImportResponseFile2='C:\Users\2013-FWQ\Desktop\AI1-01_20190513084355.txt';
% ImportResponseFile3='C:\Users\2013-FWQ\Desktop\AI1-01_20190513084355.txt';
% 
% %OutPutSourceAccTh='C:\Users\Administrator\Desktop\zhong\��3ͨ��InputAcceleration';
% OutPutResponseAccTh='C:\Users\2013-FWQ\Desktop\��1ͨ��Acceleration';
% OutPutResponsePSD='C:\Users\2013-FWQ\Desktop\��1ͨ��PSD';
% %OutPutResponseVelTh='C:\Users\Bin Peng\Desktop\201901����������\191221\��3ͨ��Velocity';
% %OutPutResponseDisTh='C:\Users\Bin Peng\Desktop\201901����������\191221\��3ͨ��Displacement';
% % OutPutXCorr='C:\Users\Bin Peng\OneDrive\Temp\��4ͨ��XCorrelation'
% %OutPutCoheret='C:\Users\Bin Peng\OneDrive\Temp\��4ͨ��Coherent';
% 
% SigLen=12803893-12;
% NumSeg=1;
% SamFreq=50;
% %ʩ������ʦ�豸
% % SenseOfAcc=(1/49)*980 ;%Gal/v
% % SenseOfInstru=23.8; %Gal/v
% % GainOfInstru=1;
% %�����豸
% SenseOfAcc=32.8 ;%Gal/v
% SenseOfInstru=34.2; %Gal/v
% GainOfInstru=1;
% R=1/((SenseOfInstru/SenseOfAcc)*GainOfInstru);
% th=[1/SamFreq:1/SamFreq:SigLen*(1/SamFreq)]';
% DtrOdr=1;
% 
% %Input data
% % importdata(ImportSourceFile1,' ',1);
% % A1=ans.data*R;
% % importdata(ImportSourceFile2,' ',1);
% % A2=ans.data*R;
% % importdata(ImportSourceFile3,' ',1);
% % A3=ans.data*R;
% % Asum=[A1;A2;A3];
% % for i=1:NumSeg;
% %    for j=1:SigLen;
% %     A(j,i)=Asum((i-1)*SigLen+j);
% %    end
% % end
% % for k=1:NumSeg;
% %      A(:,k)=detrend(A(:,k));ȥ����������
% %      A(:,k)=A(:,k)-polyval(polyfit(th,A(:,k),DtrOdr),th);%ȥһ��������
% % end
% % In_th=mean(A,2);
% 
% importdata(ImportResponseFile1,'\t',12);
% B1=ans.data*R;
% importdata(ImportResponseFile2,'\t',12);
% B2=ans.data*R;
% importdata(ImportResponseFile3,'\t',12);
% B3=ans.data*R;
% Bsum=[B1;B2;B3];
% for i=1:NumSeg
%    for j=1:SigLen
%      B(j,i)=Bsum((i-1)*SigLen+j);
%    end
% end
% for k=1:NumSeg
%      %B(:,k)=detrend(B(:,k));%��ȥ����������
%      B(:,k)=B(:,k)-polyval(polyfit(th,B(:,k),DtrOdr),th);%ȥһ��������
% end
% Acc_th=mean(B,2);
% load('C:\Users\Bin Peng\Desktop\MAT20200420-1254���ǰ�����\MAT20200420-1254ԭʼ����\202004201254.mat');
  load('C:\Users\Bin Peng\Desktop\1�Ų��\1�Ų��_2020_04_20_145953.mat');
 for i=1:24
 %   SigLen=size(Datas,1);
%   Acc_th=Datas(SigLen+1:2*SigLen,i);
%   Acc_th=detrend(Acc_th);
   SamFreq=500;
%   th=[1/SamFreq:1/SamFreq:SigLen*(1/SamFreq)]';
   OutPutResponsePSD=['C:\Users\Bin Peng\Desktop\',num2str(i)];
%% show the In_ch
% figure(1);
% th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% plot(th,In_th);
% title('Input accelartion');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Acceleration (cm/s^2)');
% ylim([1.1*min(In_th),1.1*max(In_th)]);
% print('-f1',OutPutSourceAccTh,'-painters','-dmeta','-r600');

%% show the Acc_ch
% figure(2);
% plot(th,Acc_th)
% title('Accelartion');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Acceleration (cm/s^2)');
% ylim([1.1*min(Acc_th),1.1*max(Acc_th)]);
% print('-f2',OutPutResponseAccTh,'-painters','-dmeta','-r600');

%% PSD maker of Acc_th
% figure(3);
% % Acc_f=fft(Acc_th,SigLen);
% % Psd=Acc_f.*conj(Acc_f)/SigLen;
% % f=SamFreq/SigLen*(1:(0.5*SigLen))';
% % plot(f,Psd(1:(0.5*SigLen)));
% % [Psd,f]=periodogram(Acc_th,[],[],SamFreq);
% % plot(f,Psd);
% [Psd,f]=pwelch(Acc_th,[],[],[],SamFreq);%Ĭ�ϴ����ֶ������ص�����
% plot(f,Psd);
% title('Power spectral density');
% xlabel('Frequency (Hz)');
% xlim([1,30]);
% % ylim([0,0.05]);
% %(���ṹ����ѧ��p.297�ϵĵ�λΪ|Acc|^{2}/2Hz,��˴���ͬ��The units of the PSD estimate are in squared magnitude units of the time series data per unit frequency.
% %For example, if the input data is in volts, the PSD estimate is in units of squared volts per unit frequency. 
% %For a time series in volts, if you assume a resistance of 1 �� and specify the sample rate in hertz, the PSD estimate is in watts per hertz.
% ylabel('|Acc|^{2}/Hz');
% % ylim([0,1.1*max(Psd)]);
% % [~,S] = max(Psd);
% % Smax=S*(SamFreq/SigLen);
% % sSmax=num2str(Smax,2);
% % text(1,0.9*max(Psd),['\rightarrow',sSmax]);
% print('-f3',OutPutResponsePSD,'-painters','-djpeg','-r600');
 
 
% %% Velocity time history maker of Acc_th
% figure(4);
% % th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% % th=th';
% V_th=cumtrapz(th,Acc_th)*(1/SamFreq);
% % V_th(1)=0;
% % for i=2:SigLen
% %     V_th(i)=V_th(i-1)+(Acc_th(i)+Acc_th(i-1))*(1/SamFreq)*0.5;
% % end
% [pksv,ls]=max(abs(V_th));
% plot(th,V_th);
% title('Velocity');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Velocity (cm/s)');
% ylim([1.1*min(V_th),1.1*max(V_th)]);
% Vmax=max(abs(V_th));
% sVmax=num2str(Vmax,2);
% text(ls*(1/SamFreq),max(abs(V_th))*sign(V_th(ls)),['\rightarrow',sVmax]);
% print('-f4',OutPutResponseVelTh,'-painters','-djpeg','-r600');

%% Displacement time history maker of Acc_th
% figure(5);
% th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% th=th';
% D_th=cumtrapz(th,V_th)*(1/SamFreq);
% % D_th(1)=0;
% % for j=2:SigLen
% %     D_th(j)=D_th(j-1)+(V_th(j)+V_th(j-1))*(1/SamFreq)*0.5;
% % end   
% plot(th,D_th);
% title('Displacement');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Displacement (cm)');
% ylim([1.1*min(D_th),1.1*max(D_th)]);

%% divide the In_ch to n sections and make the averaged PSD
% figure(1);
dD=reshape(Datas(:,i),300000,[]);
% DD_th=cumtrapz(th,Datas(:,i))*(1/SamFreq);
% 
for k=1:size(dD,2)
dD(:,k)=dD(:,k)-mean(dD(:,k));
dD(:,k)=detrend(dD(:,k));
end

SigLen=size(dD,1);
th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% th=th';
% D_th=cumtrapz(th,Datas(:,i))*(1/SamFreq);
% DD_th=reshape(D_th,300000,[]);

for j=1:(size(dD,2)-1)
DD_th(:,j)=cumtrapz(th,dD(:,j))*(1/SamFreq);
%[Psd(:,j),f]=pwelch(DD_th(:,j),[],[],[],SamFreq);%Ĭ�ϴ����ֶ������ص�����
%[Psd(:,j),f]=pwelch(dD(:,j),[],[],[],SamFreq);%Ĭ�ϴ����ֶ������ص�����
end

% MPsd=mean(Psd,2);
plot(th,DD_th(:,j));
title('Averaged power spectral density');
xlabel('time (s)');
xlim([0,600]);
ylabel('Dis');

% th=[0:1/SamFreq:(SigLen-1)/SamFreq];
% plot(th,In_th);
% title('Input accelartion');
% xlabel('Time (s)');
% xlim([0,SigLen/SamFreq]);
% ylabel('Acceleration (cm/s^2)');
% ylim([1.1*min(In_th),1.1*max(In_th)]);
print('-f1', OutPutResponsePSD,'-painters','-djpeg','-r600');
end
%% Cross correlation  
% figure(6);
% [acor,lag] = xcorr(In_th,Acc_th);
% acor=acor/max(acor);
% [~,tI] = max(abs(acor));
% tlagDiff = lag(tI);
% timeDiff = tlagDiff/SamFreq;
% plot(lag/SamFreq,acor);
% title('Cross correlation function');
% xlim([-0.5*((length(lag)-1)/SamFreq),0.5*((length(lag)-1)/SamFreq)]);
% xlabel('Time lag (s)');
% ylim([1.1*min(acor),1.1*max(acor)]);
% ylabel('|Acc|^2');

%% Cross coherent  
% figure(7);
% fs=[0:1/SamFreq:SamFreq/2];
% [Cxy,F] = mscohere(In_th,Acc_th,[],[],fs,SamFreq);
% [Cxy,F] = mscohere(In_th,Acc_th,[],[],[],SamFreq);%Ĭ�ϴ����ֶ������ص�����
% [~,fI] = max(abs(Cxy));
% flagDiff = F(fI);
% frequencyDiff = flagDiff/SamFreq;
% plot(F,Cxy);
% title('Cross coherent function');
% xlim([0,(length(F)-1)/SamFreq]);
% xlabel('frequency (Hz)');
% ylim([0,1]);
% print('-f7',OutPutCoheret,'-painters','-dmeta','-r600');

