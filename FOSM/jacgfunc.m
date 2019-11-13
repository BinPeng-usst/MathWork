function g_out = jacgfunc( jacg , x_in ) %功能函数偏导数计算即雅可比矩阵计算 
f = x_in(1); 
c = x_in(2); 
V = x_in(3); 
for i = 1:3 
   g_out(i) = eval(jacg(i));%1为对x的导数,2为对y的导数 
end
%将以上内容保存为jacgfunc.m 
