clear all;
for i=1:3
  I=num2str(i);
  for j=1:7
  J=num2str(j);
  str=['G:\������Ŀ\�����°뵼�幫˾��������\180130\8#\',I,'��ͨ�� ',J,'.dat'];
  temp=dlmread(str,' ',1,0);
  eval(['var',I,J,'=','temp',';']);
  %load(str);
  end
end

for m=1:3
    M=num2str(m);
    for n=1:7
        N=num2str(n);
        eval(['P','=','var',M,N,';']);
        fig=figure('Name',['var',M,N]);
        plot(P);
        print(fig,['H:\MathWork\matlab\','var-',M,N],'-painters','-dmeta','-r600');
    end 
end