function [X,Y,Z] = TubeLike(x,y,z,r)
% 绘制三维管道型立体
% TubeLike(x,y,z)    绘制三维管道型几何体，输入参数x,y,z分别为管道中心线各点处
%                    的坐标，x,y和z应为等长的向量，此时管道具有统一的半径1.
%
% TubeLike(x,y,z,r)  输入参数r用来指定管道半径。r可以是标量，也可以是与x,y,z等
%                    长的向量。当r是标量时，管道具有统一半径；当r是向量时，管道
%                    各截面处可以有不同的半径。
%
% TubeLike(x)        输入参数x为3行或3列的矩阵，用来指定管道中心线各点处的坐标，
%                    此时管道具有统一的半径1.
%
% TubeLike(x,r)      输入参数x为3行或3列的矩阵，用来指定管道中心线各点处的坐标，
%                    输入参数r（标量或向量）用来指定管道半径。
%
% [X,Y,Z] = TubeLike(...)  输出三维管道型几何体的网格数据X，Y和Z.
%                   
% CopyRight:xiezhh（谢中华）
% 2011.6.26
%
% Example:
%     t = linspace(0,2*pi,50);
%     x = sin(t);
%     y = cos(t);
%     z = cos(t/2); 
%     r = sin(t);
%     TubeLike(x,y,z,r)


if nargin>=1 && nargin<=2
    % 检查数据维数是否正确
    [m,n] = size(x);
    p = min(m,n);  % 维数
    if p ~= 3
        error('应输入三维样本数据,并且样本容量应大于3');
    end
    % 把样本观测值矩阵转置，使得行对应变量，列对应观测
    if m >= n
        x = x';
    end 
    yd = x(2,:);
    zd = x(3,:);
    xd = x(1,:);
    if nargin == 1
        r = ones(size(xd));
    else
        if isvector(y)
            if numel(y) == 1
                r = y*ones(size(xd));
            elseif numel(y)>1 && numel(y) == length(xd)
                r = y(:)';
            else
                error('半径应为标量或与x等长的向量');
            end
        else
            error('半径应为标量或与x等长的向量');
        end
    end
elseif nargin>=3 && nargin<=4
    if isvector(x) && isvector(y) && isvector(z)
        numxyz = [numel(x),numel(y),numel(z)];
        if any(numxyz-min(numxyz))
            error('管道中心坐标x,y,z应为等长的向量');
        else
            xd = x(:)';
            yd = y(:)';
            zd = z(:)';
        end
    else
        error('管道中心坐标x,y,z应为等长的向量');
    end
    if nargin == 3
        r = ones(size(xd));
    else
        if isvector(r)
            if numel(r) == 1
                r = r*ones(size(xd));
            elseif numel(r)>1 && numel(r) == length(xd)
                r = r(:)';
            else
                error('半径应为标量或与x等长的向量');
            end
        else
            error('半径应为标量或与x等长的向量');
        end
    end
else
    error('至少需要1个输入参数，至多需要4个输入参数');
end

t = linspace(0,2*pi,30)';  % 角度向量

% x的一阶差分
dx = diff(xd);
dx = [dx(end) dx]; 
% y的一阶差分
dy = diff(yd);
dy = [dy(end) dy]; 
% z的一阶差分
dz = diff(zd);
dz = [dz(end) dz]; 

%计算法线与y轴正向夹角余弦，法线与y轴正向夹角正弦负值
den1 = sqrt(dx.^2 + dy.^2);
cy = dy./den1;
cy(den1 == 0) = 1;
sy = -dx./den1;
sy(den1 == 0) = 0;

% 计算法线与z轴正向夹角余弦，法线与z轴正向夹角正弦负值
den2 = sqrt(dx.^2 + dy.^2 + dz.^2);
cz = dz./den2;
cz(den2 == 0) = 1;
sz = -sqrt(dx.^2 + dy.^2)./den2;
sz(den2 == 0) = 0;

OneMat = ones(numel(t),1);  % 1向量

% 管道中心线坐标矩阵
Xcenter = OneMat*xd;
Ycenter = OneMat*yd;
Zcenter = OneMat*zd;

% 单位圆坐标数据
x0 = cos(t);
y0 = sin(t);

% 计算三维管道型几何体的网格数据X，Y和Z
Xgrid = Xcenter + x0*(r.*cy) - y0*(r.*cz.*sy);
Ygrid = Ycenter + x0*(r.*sy) + y0*(r.*cz.*cy);
Zgrid = Zcenter + y0*(r.*sz);

% 输出图形或网格数据
if nargout == 0
    surf(Xgrid,Ygrid,Zgrid);
else
    X = Xgrid;
    Y = Ygrid;
    Z = Zgrid;
end

