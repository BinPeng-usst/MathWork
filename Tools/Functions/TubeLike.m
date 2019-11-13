function [X,Y,Z] = TubeLike(x,y,z,r)
% ������ά�ܵ�������
% TubeLike(x,y,z)    ������ά�ܵ��ͼ����壬�������x,y,z�ֱ�Ϊ�ܵ������߸��㴦
%                    �����꣬x,y��zӦΪ�ȳ�����������ʱ�ܵ�����ͳһ�İ뾶1.
%
% TubeLike(x,y,z,r)  �������r����ָ���ܵ��뾶��r�����Ǳ�����Ҳ��������x,y,z��
%                    ������������r�Ǳ���ʱ���ܵ�����ͳһ�뾶����r������ʱ���ܵ�
%                    �����洦�����в�ͬ�İ뾶��
%
% TubeLike(x)        �������xΪ3�л�3�еľ�������ָ���ܵ������߸��㴦�����꣬
%                    ��ʱ�ܵ�����ͳһ�İ뾶1.
%
% TubeLike(x,r)      �������xΪ3�л�3�еľ�������ָ���ܵ������߸��㴦�����꣬
%                    �������r������������������ָ���ܵ��뾶��
%
% [X,Y,Z] = TubeLike(...)  �����ά�ܵ��ͼ��������������X��Y��Z.
%                   
% CopyRight:xiezhh��л�л���
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
    % �������ά���Ƿ���ȷ
    [m,n] = size(x);
    p = min(m,n);  % ά��
    if p ~= 3
        error('Ӧ������ά��������,������������Ӧ����3');
    end
    % �������۲�ֵ����ת�ã�ʹ���ж�Ӧ�������ж�Ӧ�۲�
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
                error('�뾶ӦΪ��������x�ȳ�������');
            end
        else
            error('�뾶ӦΪ��������x�ȳ�������');
        end
    end
elseif nargin>=3 && nargin<=4
    if isvector(x) && isvector(y) && isvector(z)
        numxyz = [numel(x),numel(y),numel(z)];
        if any(numxyz-min(numxyz))
            error('�ܵ���������x,y,zӦΪ�ȳ�������');
        else
            xd = x(:)';
            yd = y(:)';
            zd = z(:)';
        end
    else
        error('�ܵ���������x,y,zӦΪ�ȳ�������');
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
                error('�뾶ӦΪ��������x�ȳ�������');
            end
        else
            error('�뾶ӦΪ��������x�ȳ�������');
        end
    end
else
    error('������Ҫ1�����������������Ҫ4���������');
end

t = linspace(0,2*pi,30)';  % �Ƕ�����

% x��һ�ײ��
dx = diff(xd);
dx = [dx(end) dx]; 
% y��һ�ײ��
dy = diff(yd);
dy = [dy(end) dy]; 
% z��һ�ײ��
dz = diff(zd);
dz = [dz(end) dz]; 

%���㷨����y������н����ң�������y������н����Ҹ�ֵ
den1 = sqrt(dx.^2 + dy.^2);
cy = dy./den1;
cy(den1 == 0) = 1;
sy = -dx./den1;
sy(den1 == 0) = 0;

% ���㷨����z������н����ң�������z������н����Ҹ�ֵ
den2 = sqrt(dx.^2 + dy.^2 + dz.^2);
cz = dz./den2;
cz(den2 == 0) = 1;
sz = -sqrt(dx.^2 + dy.^2)./den2;
sz(den2 == 0) = 0;

OneMat = ones(numel(t),1);  % 1����

% �ܵ��������������
Xcenter = OneMat*xd;
Ycenter = OneMat*yd;
Zcenter = OneMat*zd;

% ��λԲ��������
x0 = cos(t);
y0 = sin(t);

% ������ά�ܵ��ͼ��������������X��Y��Z
Xgrid = Xcenter + x0*(r.*cy) - y0*(r.*cz.*sy);
Ygrid = Ycenter + x0*(r.*sy) + y0*(r.*cz.*cy);
Zgrid = Zcenter + y0*(r.*sz);

% ���ͼ�λ���������
if nargout == 0
    surf(Xgrid,Ygrid,Zgrid);
else
    X = Xgrid;
    Y = Ygrid;
    Z = Zgrid;
end

