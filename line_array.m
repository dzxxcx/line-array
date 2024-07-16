%% 一维线阵归一化方向图
% 归一化方向图+目标最大波束指向

clc
clear all
close all
eps = 0.00001;
N = 8;                 % 阵元数量
f = ones(1,N);         % 阵元的方向图函数
I = ones(1,N);         % 激励电流
d_lambda = 1/2;        % 阵元间隔/波长
thetab = 30;           % 目标最大波束指向(°)
deltaphib = 2*pi*d_lambda*sind(thetab);      
% deltaphib = 0*pi/180 % 激励电流控制的初相(°) -> (rad)
theta = -80:1:80;      % 俯仰角角度
phi =  -180:1:180;     % 方位角角度

%% 三维图
F = zeros(length(theta),length(phi));
for i = 0:N-1
    Fi= f(i+1)*I(i+1)*exp(-j*i*deltaphib)*exp(j*i*2*pi*d_lambda*cosd(theta)'*sind(phi));
    F = F + Fi;
end
F = F / N; % 归一化方向图
subplot(2,2,1);
mesh(phi,theta,abs(F));
xlabel('方位角/\circ');ylabel('俯仰角/\circ');zlabel('归一化方向图');
title('三维图');

%% 三维截面图
theta0 = zeros(1,length(theta));      % 俯仰角角度
F0 = zeros(length(theta0),length(phi));
for i = 0:N-1
    Fi0= I(i+1)*exp(-j*i*deltaphib)*exp(j*i*2*pi*d_lambda*cosd(theta0)'*sind(phi));
    F0 = F0 + Fi0;
end
F0 = F0 / N; % 归一化方向图
subplot(2,2,2);
mesh(phi,theta0,abs(F0));
xlabel('方位角/\circ');ylabel('俯仰角/\circ');zlabel('归一化方向图');
title('三维截面图');

%% 二维截面图
subplot(2,2,3);
plot(phi,abs(F0(1,:)));
xlabel('方位角/\circ');ylabel('归一化方向图');
title('二维截面图');
ylim([0, max(abs(F0(1,:)))]);
grid;

%% 极坐标图
subplot(2,2,4);
polarplot(phi*pi/180,abs(F(floor(length(theta)/2)+1,:)));
title('极坐标图');
