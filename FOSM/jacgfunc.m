function g_out = jacgfunc( jacg , x_in ) %���ܺ���ƫ�������㪤���ſɱȾ������ 
f = x_in(1); 
c = x_in(2); 
V = x_in(3); 
for i = 1:3 
   g_out(i) = eval(jacg(i));%1Ϊ��x�ĵ���,2Ϊ��y�ĵ��� 
end
%���������ݱ���Ϊjacgfunc.m 
