%{ 
**************************************************************************************************************************************************
**************************************************************************************************************************************************
【起始时间】：2023.6.19
【参与人】：杨浩东
【参考文献】：
【中文】
[1]王肖,郭杰,唐胜景,祁帅.基于解析剖面的时间协同再入制导[J].航空学报,2019,40(03):239-250.
[2]姜鹏,郭栋,韩亮,李清东,任章.多飞行器再入段时间协同弹道规划方法[J].航空学报,2020,41(S1):171-183.
[3]乔浩,李师尧,李新国.多高超声速飞行器静态协同再入制导方法[J].宇航学报,2020,41(05):541-552.
[4]许强强. 多飞行器协同轨迹规划与制导研究[D].国防科技大学,2020.DOI:10.27052/d.cnki.gzjgu.2020.000380.
[5]许强强. 高超声速多飞行器编队协同动力学研究[D].国防科学技术大学,2016.
[6]乔浩. 高超声速滑翔导弹协同制导技术研究[D].西北工业大学,2019.DOI:10.27406/d.cnki.gxbgu.2019.000328.
【英文】
[1]Jianglong Y U, Xiwang D, Qingdong L I, et al. Cooperative guidance strategy for multiple hypersonic gliding vehicles system[J]. Chinese Journal of Aeronautics, 2020, 33(3): 990-1005.
[2]Yu W, Yao Y, Chen W. Analytical cooperative entry guidance for rendezvous and formation flight[J]. Acta astronautica, 2020, 171: 118-138.
[3]Li Z, He B, Wang M, et al. Time-coordination entry guidance for multi-hypersonic vehicles[J]. Aerospace Science and Technology, 2019, 89: 123-135.
[4]Wang H, Guo J, Wang X, et al. Time-coordination entry guidance using a range-determined strategy[J]. Aerospace Science and Technology, 2022, 129: 107842.
[5]Liang Z, Yu J, Ren Z, et al. Trajectory planning for cooperative flight of two hypersonic entry vehicles[C]//21st AIAA International Space Planes and Hypersonics Technologies Conference. 2017: 2251.
[6]Chu H, Li J, Dong Y, et al. Improved MPSP Method-based Cooperative Re-entry Guidance for Hypersonic Gliding Vehicles[C]//Matec Web of Conferences. EDP Sciences, 2017, 114: 01002.
[7]Yu W ,  Chen W ,  Jiang Z , et al. Analytical entry guidance for coordinated flight with multiple no-fly-zone constraints[J]. Aerospace science and technology, 2019, 84(JAN.):273-290.
【说明】：
1、首先打算先参考我本科毕业论文，重新修改实现
2、基于已实现的内容再进一步寻找创新点
**************************************************************************************************************************************************
**************************************************************************************************************************************************
%}
%% 初始化设置
close all
% clearvars -except XX
clear
clc

% 获取当前文件位置，拆分，获取当前文件所在路径，将所在路径及其子文件加入路径中
% % 将当前级的上一级整体作为路径
% currentFolder = pwd;
% [filepath,name,ext] = fileparts(currentFolder);
% addpath(genpath(filepath))

% 将当前级作为路径
currentFolder = pwd;
addpath(genpath(currentFolder))

%% 禁用奇异矩阵警告，不然太多命令行窗口信息了
warning('off')
% warning('off','MATLAB:singularMatrix')
% warning

%% 初始化
Parameter_CAV_H
%% 输入
% 起终点速度、高度
V0 = 7000; H0 = 54e3;
Vf = 2000; Hf = 25e3;
HV0f = [H0,Hf;V0,Vf];

r0 = H0 + c.Re;
lambda0 = 0;
phi0 = 0;
theta0 = 0;
psi0 = 45*pi/180;
S0 = 0;
X = [r0,lambda0,phi0,V0,theta0,psi0,S0];

% 攻角用角度制，因为在气动模型中直接带入的角度制
alpha = 20;
beta = 0;
sigma = 10*pi/180;
u = [alpha,beta,sigma];

Target = zeros(2,1);
% lambdaT = 55;
% phiT = 55;
lambdaT = 50;
phiT = 50;
Target(1) = lambdaT*pi/180;
Target(2) = phiT*pi/180;

%% 基于大圆弧假设的初始预期航程计算
latT = phiT; lonT = lambdaT;
distT_deg = distance(phi0,lambda0,latT,lonT);
% distT_s = distT_deg*pi/180*c.Re + 2.6069e+05;
% distT_s = distT_deg*pi/180*c.Re*(1+1/75); % 50经度 50纬度 10航向角走廊
distT_s = distT_deg*pi/180*c.Re;
% distT_s = sqrt((distT_deg*pi/180*c.Re)^2 + (H0 - Hf)^2);
%% 基于给定时间、航程的高度系数求解
% 参考，1200s基本对应6500km
% 之后给定航程可设置的时间参考这个标准上下浮动
t0 = 1550;
% t0 = 1600;
% s0 = 6.5e6;
s0 = distT_s;
K0 = [0.5,0.5];

A = [];
b = [];
Aeq = [];
beq = [];
% lb = [0.2,0.2];
% ub = [0.8,0.8];
lb = [0.2,0];
ub = [1,1];
nonlcon = [];
% options = optimoptions("fmincon","Algorithm","sqp", ...
%     'TolFun', 1e-5,'TolX', 1e-5,'TolCon',5e-5);
options = optimoptions("fmincon","Algorithm","sqp");


V = V0; % 当前时刻速度
H = H0;
[K,f] = fmincon(@(K) fun_Object(K,t0,s0,HV0f,c,V,H),K0,A,b,Aeq,beq,lb,ub,nonlcon,options);

% 控制点速度
V1 = Vf + 1*(V0 - Vf)/4;
V2 = Vf + 2*(V0 - Vf)/4;

% 控制点高度设计
% 高度系数
% K = [0.5,0.5];
k1 = K(1); 
k2 = K(2);
[H1min,H1max] = fun_FlightCorridor(c,V1);
H1 = H1min + k1*(H1max - H1min);
[H2min,H2max] = fun_FlightCorridor(c,V2);
H2 = H2min + k2*(H2max - H2min);

%% 基于H-V剖面的积分计算航程、时间
HV = [H0,H1,H2,Hf;V0,V1,V2,Vf];
TS0 = [0 0];
Vspan = [Vf V0];
[Vts_list,TS_list] = ode45(@(V,TS) fun_TSpre_old(c,HV,V),Vspan,TS0);
sgo1 = TS_list(end,2); % 基于H-V积分计算的剩余航程
tgo1 = TS_list(end,1); % 基于H-V积分计算的剩余航程

%% 标准剖面高度
[Hpro,~,~] = fun_HVprofile(HV,V);

%% List
Listmax = 1e3;
list_X = zeros(Listmax,7); 
list_u = zeros(Listmax,3);
list_t = zeros(Listmax,1);
list_K = zeros(Listmax,2);
list_f = zeros(Listmax,1);
list_Hpro = zeros(Listmax,1);
list_tk = zeros(1,1);

%% 仿真主函数
nmax = 1e5;
dt = 0.1;
t = 0;
nk = 0; % 高度系数更新相关计数器
tStart = tic;
for n = 1:nmax
    %% List记录
    list_X(n,:) = X;
    list_u(n,:) = u;
    list_t(n,1) = t;
    list_K(n,:) = K;
    list_f(n,1) = f;
    list_Hpro(n,1) = Hpro;

    %% 终止条件
    if X(4) < Vf
        break
    end
    
    %% 若干次制导循环更新一次高度系数
    nk = nk + 1;
    % if ((t <= 800) && (nk == 1000)) || ((t > 800) && (t <= 1200) && (nk == 100)) || ((t > 1200) && (t <= 1400) && (nk == 10)) || ((t > 1400) && (nk == 1))
    if nk == 1000
        tic
        %% 高度系数校正
        H = X(1) - c.Re;
        lambda = X(2)*180/pi; % 当前时刻状态
        phi = X(3)*180/pi;
        V = X(4);
        sgo_deg = distance(phi,lambda,latT,lonT);
        sgo = sgo_deg*pi/180*c.Re*(1+1/75);
        % sgo = sgo_deg*pi/180*c.Re*(2*H+c.Re)/c.Re;
        % sgo = sgo_deg*pi/180*c.Re;
        % sgo = sqrt((sgo_deg*pi/180*c.Re)^2 + (H - Hf)^2);
        tgo = t0 - t;
        [K,f] = fun_Kdetermination(c,tgo,sgo,HV0f,V,H,K);
        % if f > list_f(n)
        %     K = list_K(n,:);
        % end
        %% 基于高度系数的控制点设计
        % V1 = Vf + (V - Vf)/4;
        % V2 = Vf + 2*(V - Vf)/4;
        k1 = K(1);
        k2 = K(2);
        % [H1min,H1max] = fun_FlightCorridor(c,V1); % V1如果是固定的，则不用每次都重新算
        H1 = H1min + k1*(H1max - H1min);
        % [H2min,H2max] = fun_FlightCorridor(c,V2); % V2如果是固定的，则不用每次都重新算
        H2 = H2min + k2*(H2max - H2min);
        HV = [H0,H1,H2,Hf;V0,V1,V2,Vf];
        % HV = [H,H1,H2,Hf;V,V1,V2,Vf];

        nk = 0;
        tk = toc;
        list_tk = [list_tk;tk];
    end
    

    V = X(4);
    [Hpro,~,~] = fun_HVprofile(HV,V);

    u = fun_Guidance(c,X,u,HV,Target);
    %% RK4
    dX1 = fun_ReentryDynamics(c,X,u);
    dX2 = fun_ReentryDynamics(c,X + 0.5*dt*dX1,u);
    dX3 = fun_ReentryDynamics(c,X + 0.5*dt*dX2,u);
    dX4 = fun_ReentryDynamics(c,X + dt*dX3,u);
    X = X + dt*(dX1 + 2*dX2 + 2*dX3 + dX4)/6;
    t = t + dt;

end
tEnd = toc(tStart);

%% 基于大圆弧假设的实际航程计算
lat1 = 0; lon1 = 0;
lat2 = X(3); lon2 = X(2);

dist1_deg = distance(0,0,lat2*180/pi,lon2*180/pi);
dist1_s = dist1_deg*pi/180*c.Re;
dist2_deg = acos(sin(lat2)*sin(lat1) + cos(lat2)*cos(lat1)*cos(lon2-lon1))*180/pi;
dist2_s = dist2_deg*pi/180*c.Re;
%%
Plot_Corridor
Plot_State

%% 解析时间预测
% 比积分计算快几个数量级
L_D = 2; % 升阻比
Vcr = sqrt(c.g*c.Re); % 第一宇宙速度
v = 0; % 倾侧角

% [1]韩嘉俊,王小虎,郝昀,张后军.带有时间约束的再入滑翔轨迹设计[J].宇航学报,2020,41(04):438-446.
% 基于航程计算剩余时间，基于解析航程算的不太准，倒是先积分算航程再代入挺准的
sgo2 = c.Re/2*L_D*cos(v)*log((Vcr^2 - Vf^2)/(Vcr^2 - V0^2));
tgo22 = sgo2/Vcr*log((Vcr + V0)*(Vcr - Vf)/(Vcr - V0)/(Vcr + Vf))/log((c.g*c.Re - Vf^2)/(c.g*c.Re - V0^2));
tgo21 = sgo1/Vcr*log((Vcr + V0)*(Vcr - Vf)/(Vcr - V0)/(Vcr + Vf))/log((c.g*c.Re - Vf^2)/(c.g*c.Re - V0^2));


sgo3 = c.Re/2*L_D*log((c.g*(H0+c.Re) - Vf^2)/((c.g*(H0+c.Re) - V0^2)));

% [1]王洁瑶,江涌,钟世勇,熊灵芳.高超声速远程导弹弹道解析估算与特性分析[J].宇航学报,2016,37(04):403-410.
% 基于常升阻比计算剩余时间
tgo3 = L_D*Vcr/(2*c.g)*(log((Vcr-Vf)/(Vcr+Vf)) - log((Vcr-V0)/(Vcr+V0)));

%% 误差
% 时间误差
% clc
Terror = t0 - t;
if Terror > 0
    fprintf('终端时间少了 %4.2f s\n',abs(Terror))
else
    fprintf('终端时间多了 %4.2f s\n',abs(Terror))
end

% 航程误差
lambda = X(2)*180/pi; % 当前时刻状态
phi = X(3)*180/pi;
Serror = distance(phi,lambda,latT,lonT)*pi/180*c.Re;
if (lambda < lambdaT) && (phi < phiT)
    fprintf('终端航程少了 %4.2f km\n',Serror/1e3)
else
    fprintf('终端航程多了 %4.2f km\n',Serror/1e3)
end

% 高度误差
Herror = X(1) - c.Re - Hf;
if Herror > 0
    fprintf('终端高度多了 %4.2f m\n',Herror)
else
    fprintf('终端高度少了 %4.2f m\n',Herror)
end

% 速度误差
Verror = X(4) - Vf;
if Verror > 0
    fprintf('终端速度多了 %4.4f m/s\n',Verror)
else
    fprintf('终端速度少了 %4.4f m/s\n',Verror)
end
%%
% 真实航程
fprintf('真实航程为 %4.2f km\n',X(7)/1e3)
% 初始大圆弧预测航程
fprintf('初始大圆弧预测航程为 %4.2f km\n',distT_deg*pi/180*c.Re/1e3)
% 考虑补偿的大圆弧预测航程
fprintf('考虑补偿的大圆弧预测航程为 %4.2f km\n',distT_s/1e3)

fprintf('真实航程与预测航程之差为 %4.2f km\n',X(7)/1e3-distT_deg*pi/180*c.Re/1e3)
fprintf('真实航程与补偿航程之差为 %4.2f km\n',X(7)/1e3-distT_s/1e3)

