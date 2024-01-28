clear;clc;close all
f=@(i,j) (i-j)+i*(i==j);
% myColor = [0 0 1;
%     0 0.5 0
%     1 0 0
%     0 0.75 0.75
%     0.75 0 0.75
%     0.75 0.75 0
%     0.25 0.25 0.25];
%                             2022.10.14
% 现在的版本在无DO的情况下也能满足队列稳定                 qin的三个增益可能需要调整以适应拓扑
%  自己的控制器，与无DO对比，同时与其他控制器对比，增加了正弦扰动（效果还行）   DO参数还得再调一下
%% 过程控制
doPic = true;
doFeedforward = false; %是否前馈
isDelayed = false; %前馈中的邻域控制输入是否有时延
isRandom = 0 ; %1：随机初始状态,0：零偏差初始状态
isPseudoRandom = false; %1：固定的随机初始状态，0：非固定的随机初始状态
isHeterogeneous = true; %性能指标是否异质
doDesiredAccelerationInput = false; %领航车运行方式，true为期望加速度输入，false为实际加速度输入
doAccelerate = false; %领航车是否始终加速
IndexTopo = 1 ; % 拓扑类型 1.PF 2.PLF 3.TPF 4.TPLF 5.MF
flagParaMismatch = 0; % 是否存在参数不匹配
flagCtrlSgn2Sat = 0  ; % sat 替代 sgn
flagCtrlSat = 1 ; % 是否存在输入饱和
flagmisdis = 1 ;  %是否存在失配扰动dd1,dd2
flagDO = 1  ;   %是否存在扰动观测器
flagfiniteDO = 1 ; %   1有限时间DO,  0线性wangDO，同时切换为wang2020控制器
flagFinitefTimeCtrl = 1 ; % 采用有限时间控制，0为线性反馈控制（qin鲁棒控制器）
%% 初始化
ddes= 8 ;% 最小安全距离。。。。6
v0= 16 ; %初始车速   6
v00 = 16;%领航车初始车速 6
D0=[ddes;0;0];%状态偏移
h= 0.5 ;%%时距
d0=ddes+h*v0; %初始车距  12。。。。9
C=diag([1,1,1]);
k11 = 2.5; k22 = 1.26; k33 = 0.23 ;   %%%   线性反馈控制增益
% k11 = 0.1581; k22 = 0.5142; k33 = 0.2571 ;   %%%   线性反馈控制增益  重新计算
% k11 = 0.84; k22 = 2.41; k33 = 1.16 ;   %%%   线性反馈控制增益  commen
% k11 = 1; k22 = 0.5; k33 = 0;
k1 = 0.15 ;   % 对应论文中的k_i   0.6
epsilon = 1   ; %   0.3  0.5  这个参数论文中没有！！！！！！！！！取为1
delta = 0.5 ; %    对应论文中的gamma  0.45
q= 0.8 ;%耦合滑模参数0.7，对应论文中的betta
ppp = 0.3 ;  %% 0.4  对应论文中的alfa_iS
% D0 = 10;
g = 9.81;
Dmz = 5;
Dcz = 150;
am = -3; %minimum acceleration
aM = 3; %maximum acceleration
lanmda10=  10  ;  %%原 10   5
lanmda11=  8  ;%%%原 8     6
lanmda12=  7  ;%% 原 7     4
lanmda20= 5 ; %原  5     6    (((5   5
lanmda21= 6 ; % 原 6     8
lip1= 0.5 ;  %0.35
lip2= 1.5 ;  %% 0.6
%ppp = 0.3;   %%原 gamma2*gamma3/(2*gamma3 - gamma2)  wang2020参数，对应的是ppp
tstart= 1000; %绘图开始时间
tdd1 = 20+tstart/100  ;   %扰动1作用时刻
tdd2= 20+tstart/100   ;
fdd1= 0.3 ;
fdd2= 1.5 ;  %%%扰动幅值
tdd3 = 20+tstart/100  ;   %扰动3作用时刻
tdd4= 20+tstart/100   ;
fdd3= 0 ;
fdd4= 0 ;  %%%扰动幅值
l=0.1;  %wang DO参数

Nveh= 6 ; %车辆数
if Nveh<3
    msgbox('Nveh>=3')
end

Tstep=0.01; %仿真步长
deltaT = Tstep;
Nstep=tstart+5500; %仿真步数
tau = 0.3;


K=zeros(3,Nveh);%控制增益
G=zeros(Nveh,3,3);%离散状态矩阵
F=zeros(Nveh,3,1);%离散输入矩阵
x=zeros(Nveh,Nstep,3);%车辆状态
N=zeros(Nveh,3,2);%扰动系数
dd=zeros(Nveh,Nstep,2);%离散扰动
s=zeros(Nveh,Nstep);%离散滑模
S=zeros(Nveh,Nstep);%离散耦合滑模
x_wave_neighbor_average=zeros(Nveh,Nstep,3);%车辆状态误差
u_neighbor_average=zeros(Nveh,Nstep); %车辆输入
e=zeros(Nveh,Nstep,2);%跟踪误差
Vi_sum=zeros(Nveh,Nstep); %邻车速度之和
%以下为观测器部分参数
p_hat=zeros(Nveh,Nstep);
deltai1=zeros(Nveh,Nstep);
dd1_hat=zeros(Nveh,Nstep);
dd1_hat_dot=zeros(Nveh,Nstep);
dd11_hat=zeros(Nveh,Nstep);
dd11_hat_dot=zeros(Nveh,Nstep);
v_hat=zeros(Nveh,Nstep);
deltai2=zeros(Nveh,Nstep);
dd2_hat=zeros(Nveh,Nstep);
dd2_hat_dot=zeros(Nveh,Nstep);
p_hat(1,1)=d0*Nveh+h*v0-d0-0.06;


z=zeros(Nveh,Nstep);% wang DO 中间量  也就是论文中的z_hat
z_dot=zeros(Nveh,Nstep);

for i=2:Nveh
    p_hat(i,1)=p_hat(i-1,1)-d0;
end
v_hat(:,1)=v_hat(:,1)+v0;%观测器初始化

%% 拓扑
M=zeros(Nveh,Nveh); %邻接矩阵
P=zeros(Nveh,Nveh); %牵引矩阵

switch IndexTopo
    case 1 %PF
        M(2:Nveh+1:Nveh^2)=1;
        P(1,1)=1;
    case 2 %PLF
        M(2:Nveh+1:Nveh^2)=1;
        P=diag(ones(1,Nveh));
    case 3 %TPF
        M(2:Nveh+1:Nveh^2)=1;
        M(3:Nveh+1:Nveh^2-Nveh)=1;
        P(1,1)=1;P(2,2)=1;
    case 4 %TPLF
        M(2:Nveh+1:Nveh^2)=1;
        M(3:Nveh+1:Nveh^2-Nveh)=1;
        P=diag(ones(1,Nveh));
    case 5 %MF
        s=round(Nveh/2);
        M(1:s-1,2:s)=diag(ones(1,s-1));
        M((s-1)*Nveh+s+1:Nveh+1:Nveh^2)=1;
        P(s,s)=1;
        clear s
    otherwise
        msgbox('error');
end

MP=M+P; %图的特征矩阵
% DP=(MP(:,:)>0)*ones(Nveh,1); %D+P
DP=ones(Nveh,1);
Ai=[0 1 0; 0 0 1; 0 0 0];
Bi=[0 1 0]';
Ci=[1 0;0 1;0 0];

if isHeterogeneous
    r=1+(1:Nveh)/5;
    Q=[3+(1:Nveh)/5;2+(1:Nveh)/5;1+(1:Nveh)/5];
else
    r=ones(1,Nveh);
    Q=[3*ones(1,Nveh);2*ones(1,Nveh);ones(1,Nveh)];
end

randp = [0.3192    0.3129   -0.8649   -0.0301   -0.1649    0.6277    1.0933]; % 伪随机位置扰动
randv = [1.4384    0.3252   -0.7549    1.3703   -1.7115   -0.1022   -0.2414]; % 伪随机速度扰动
ctrlPara = [0.9631    0.0000   -0.0175
    0.8862    0.0000    0.0058
    1.0205    0.0000   -0.0081
    1.0134    0.0000    0.0195
    1.0082   -0.0000    0.0162
    0.9348   -0.0000   -0.0309
    1.0228    0.0000    0.0130]'; %
intOfUr = zeros(Nveh,1);
arvTimeOptm = zeros(Nveh,1);
arvTimeCtrl = zeros(Nveh,1);

for i=1:Nveh
    ri=r(i);
    Qi=diag(Q(:,i));
    
    %     Pi=are(Ai,Bi/ri*Bi',Qi);
    %     K(:,i)=Bi'*Pi/ri;
    K(:,i)=[k11,k22,k33]';
  %  K(:,1)=[2*k11,k22,k33]';
    
    if isRandom
        if isPseudoRandom
            x(i,1,:)=[randp(i)-d0*i; v0+randv(i) ;0];%初始状态
        else
            x(i,1,:)=[normrnd(0,1)-d0*i; v0+normrnd(0,1) ;0];%初始状态
        end
    else
        x(i,1,:)=[-d0*i; v0;0];%初始状态
    end
    
    syms t
    Gi=expm(Ai*Tstep); % G1=eye(3)+A*deltaT; % approximate
    Fi=int(expm(Ai*t),t,0,Tstep)*Bi;
    Fi=double(Fi); % H1=B*deltaT; % approximate
    Ni=int(expm(Ai*t),t,0,Tstep)*Ci;
    Ni=double(Ni);
    clear t
    
    G(i,:,:)=Gi;
    F(i,:,:)=Fi;
    N(i,:,:)=Ni;
end
x(:,1,1) = x(:,1,1)+(d0*Nveh+h*v0);
x0=zeros(3,Nstep); %领航车状态
x0(:,1)=[d0*Nveh+h*v0;v00;0];%领航车状态初始化（第一列）
u0=zeros(1,Nstep);
u=zeros(Nveh,Nstep);

%% 闭环动力学
for n=2:Nstep %步数
    x0(:,n)=x0(:,n-1);%领航车状态更新
    
    %以下为观测器更新
    p_hat(:,n)=p_hat(:,n-1);
    deltai1(:,n)=deltai1(:,n-1);
    dd1_hat(:,n)=dd1_hat(:,n-1);
    dd1_hat_dot(:,n)=dd1_hat_dot(:,n-1);
    dd11_hat(:,n)=dd11_hat(:,n-1);
    dd11_hat_dot(:,n)=dd11_hat_dot(:,n-1);
    v_hat(:,n)=v_hat(:,n-1);
    deltai2(:,n)=deltai2(:,n-1);
    dd2_hat(:,n)=dd2_hat(:,n-1);
    dd2_hat_dot(:,n)=dd2_hat_dot(:,n-1);
    if  n<=4/Tstep+tstart   %  原  1
        x0(3,n)=0;
        u0(n)=0;
    else
        if doAccelerate %恒定加速
            x0(3,n)=1;
            u0(n)=1;
        else %匀速、减速
            if n<=9/Tstep+tstart
                x0(3,n)= 1.5 ; % 1.5
                u0(n)= 1.5 ; % 1.5
                
            elseif n<=30/Tstep+tstart
                x0(3,n)=0;
                u0(n)=0;
            elseif n<=35/Tstep+tstart
                x0(3,n)=-1;
                u0(n)=-1;
            elseif n<=60/Tstep+tstart
                x0(3,n)=0;
                u0(n)=0;
            end
        end
    end
    if doDesiredAccelerationInput %领航车采用期望加速度输入
        x0(:,n)=Gi*x0(:,n-1)+Fi*u0(n);%领航车状态方程
    else %领航车采用实际加速度输入
        x0(2,n)=x0(2,n)+x0(3,n)*Tstep;%领航车速度更新
        x0(1,n)=x0(1,n)+x0(2,n)*Tstep+1/2*x0(3,n)*Tstep^2;%领航车位置更新
    end
    if arvTimeOptm(1) == 0 && x0(1,n) >= -Dmz
        arvTimeOptm(1) = n*Tstep;
        arvTimeOptm = arvTimeOptm(1) + (1:Nveh)'*d0/v00;
    end
    
    for i=1:Nveh %MP矩阵的行,即第n辆车
        Xi=reshape(x(i,n-1,:),3,1);%自车状态
        Xi_neighbor=zeros(3,1); %邻车状态
        Xi_error_sum=zeros(3,1); %状态误差
        Ui_neighbor=0; %邻车输入
        Ui_sum=0; %邻车输入之和
        Vi_neighbor=0; %邻车速度
        
        if flagmisdis
            dd(i,n,1)=  dd(i,n,1)+fdd1*sin(1*n*Tstep)*exp(-(n*Tstep-tdd1)^2) +fdd3*sin(0.5*n*Tstep) ;%%  dd(i,n,1)+0.25*sin(n*Tstep)
            dd(i,n,2)= dd(i,n,2)+fdd2*sin(3*n*Tstep)*exp(-(n*Tstep-tdd2)^2)+fdd4*sin(0.5*n*Tstep) ;%扰动更新   dd(i,n,2)+0.2*sin(n*Tstep)   -0.2*i
        else
            dd(i,n,1)=  0  ;  %%  dd(i,n,1)+0.25*sin(n*Tstep)
            dd(i,n,2)= 0 ;
        end
        
        for j=1:i %MP矩阵的列,即第n个邻车  原 j=1:Nveh
            if MP(i,j)~=0 %存在邻车
                if i==j %邻车为领航车
                    Xi_neighbor=x0(:,n);
                    Xi_neighbor(1)=Xi_neighbor(1)-h*(x(i,n-1,2)+dd(i,n,1))-Tstep*v0;%%  %  -Tstep*v0是为了消除离散化导致的初始误差
                    Xi_neighbor(2)=Xi_neighbor(2);%跟踪误差已更改
                    Vi_neighbor=x0(2,n);
                    if doDesiredAccelerationInput %领航车采用期望加速度输入
                        Ui_neighbor=u0(n);
                        Ui_neighbor_delayed=u0(n-1);
                    else %领航车采用实际加速度输入
                        Ui_neighbor=x0(3,n);
                        Ui_neighbor_delayed=x0(3,n-1);
                    end
                else %邻车非领航车
                    Xi_neighbor=reshape(x(j,n-1,:),3,1);
                    %                   if i==1
                    %                      Xi_neighbor(1)=0;
                    %                   else
                    Xi_neighbor(1)=Xi_neighbor(1)-h*(x(i,n-1,2)+dd(i,n,1));
                    %                   end
                    Xi_neighbor(2)=Xi_neighbor(2) ;%跟踪误差已更改
                    Ui_neighbor=u(j,n)+dd1_hat_dot(j,n);%已根据控制器更改
                    Ui_neighbor_delayed=u(j,n-1);
                    Vi_neighbor= x(j,n,2)+dd1_hat(j,n);  %
                end
                
                Xi_error_sum = Xi_error_sum + Xi-Xi_neighbor+D0*f(i,j);
                Vi_sum(i,n)=Vi_sum(i,n) + Vi_neighbor;
                if isDelayed
                    Ui_sum=Ui_sum + Ui_neighbor_delayed;
                else
                    Ui_sum=Ui_sum + Ui_neighbor;
                end
            else %不存在邻车
                continue;
            end
        end
        Ki=K(:,i);
        e(i,n,1) = Xi_error_sum(1);
        e(i,n,2) = Xi_error_sum(2);
        
        %加入观测器
        if flagDO
            if  flagfiniteDO   %有限do
                p_hat(i,n)=p_hat(i,n-1)+x(i,n-1,2)*Tstep+deltai1(i,n-1)*Tstep;
                deltai1(i,n)=-lanmda10*lip1^(1/3)*sign(p_hat(i,n)-x(i,n-1,1))*abs(p_hat(i,n)-x(i,n-1,1))^(2/3)+dd1_hat(i,n);
                dd1_hat(i,n)=dd1_hat(i,n)+dd1_hat_dot(i,n)*Tstep;
                dd1_hat_dot(i,n)=-lanmda11*lip1^(1/2)*sign(dd1_hat(i,n)-deltai1(i,n))*abs(dd1_hat(i,n)-deltai1(i,n))^(1/2)+dd11_hat(i,n);
                dd11_hat(i,n)=dd11_hat(i,n)+dd11_hat_dot(i,n)*Tstep;
                dd11_hat_dot(i,n)=-lanmda12*lip1*sign(dd11_hat(i,n)-dd1_hat_dot(i,n));
                v_hat(i,n)=v_hat(i,n-1)+u(i,n-1)*Tstep+deltai2(i,n-1)*Tstep;
                deltai2(i,n)=-lanmda20*lip2^(1/3)*sign(v_hat(i,n-1)-x(i,n-1,2))*abs(v_hat(i,n-1)-x(i,n-1,2))^(2/3)+dd2_hat(i,n);
                dd2_hat(i,n)=dd2_hat(i,n)+dd2_hat_dot(i,n)*Tstep;
                dd2_hat_dot(i,n)=-lanmda21*lip2*sign(dd2_hat(i,n)-deltai2(i,n));
            else  %线性do
                dd2_hat(i,n)=z(i,n)-l*e(i,n,1);
                z(i,n)=z(i,n)+z_dot(i,n)*Tstep;
                z_dot(i,n)=l*(e(i,n,2)+q*(u(i,n-1)+dd2_hat(i,n)));
                p_hat(i,n)=0;
                deltai1(i,n)=0;
                dd1_hat(i,n)=0;
                dd1_hat_dot(i,n)=0;
                dd11_hat(i,n)=0;
                dd11_hat_dot(i,n)=0;
                v_hat(i,n)=0;
                deltai2(i,n)=0;
                dd2_hat(i,n)=0;
                dd2_hat_dot(i,n)=0;
                
                
            end
        else
            p_hat(i,n)=0;
            deltai1(i,n)=0;
            dd1_hat(i,n)=0;
            dd1_hat_dot(i,n)=0;
            dd11_hat(i,n)=0;
            dd11_hat_dot(i,n)=0;
            v_hat(i,n)=0;
            deltai2(i,n)=0;
            dd2_hat(i,n)=0;
            dd2_hat_dot(i,n)=0;
        end
        if doFeedforward
            Ui = Ui_sum/DP(i)-Ki'*C*Xi_error_sum/DP(i);
        else
            if flagFinitefTimeCtrl%有限时间控制
                
                sumOfNghbAcc = Ui_sum;
                ur =  k1 * e(i,n,1)  ;   %修改   线性收敛率！！！
                %                 if intOfUr(i) == 0
                %                     intOfUr(i) = 0 ;%%e(i,n,2)/epsilon
                %                 else
                intOfUr(i) = intOfUr(i) + ur * deltaT;%ur的积分
                %                 end
                s(i,n) = e(i,n,1) + epsilon * intOfUr(i);
                if i<Nveh
                    S(i,n) = q*s(i,n)-s(i+1,n-1);
                else
                    S(i,n) = q*s(i,n);
                end
                if  flagfiniteDO   %有限do  （这里不止切换DO，还切换控制器跟wang2020）
                    if flagCtrlSgn2Sat
                        Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))+ epsilon * ur - delta/q * sign(S(i,n)) * abs(S(i,n))^ppp ) / (DP(i)*h)-dd2_hat(i,n);
                    else
                        if i<Nveh
                            if n<30        Vi_sum(i+1,n-1)=x(i+1,n-1,2)+dd1_hat(i+1,n); else  ;  end
                            Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))- epsilon * ur - delta/q * sign(S(i,n)) * abs(S(i,n))^ppp ) / (DP(i)*h)-dd2_hat(i,n)-dd11_hat(i,n) ...
                                +(DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))+k1 * e(i+1,n-1,1))/(DP(i+1)*h*q);% 去掉了e_i+1的导数那一项  DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))+
                            %+(DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))- k1 * sign(e(i+1,n-1,1)) * abs(e(i+1,n-1,1))^gamma1)/(DP(i+1)*h*q);
                        else
                            Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))- epsilon * ur - delta/q * sign(S(i,n)) * abs(S(i,n))^ppp ) / (DP(i)*h)-dd2_hat(i,n)-dd11_hat(i,n);
                        end
                        
                        if flagCtrlSat
                            Ui = max(min(Ui, aM), am);
                        end
                    end
                else   %%  线性do  wang2019控制器
                    if flagCtrlSgn2Sat
                        Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))- epsilon * ur - delta/q * sign(S(i,n)) * abs(S(i,n))^ppp ) / (DP(i)*h)-dd2_hat(i,n);
                    else
                        if i<Nveh
                            Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))- epsilon * ur - delta/q * (sign(S(i,n))+S(i,n) )) / (DP(i)*h)-dd2_hat(i,n)-dd11_hat(i,n) ...
                                +(DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))+k1 * sign(e(i+1,n-1,1)) * abs(e(i+1,n-1,1))^ppp)/(DP(i+1)*h*q);%
                            %+(DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))- k1 * sign(e(i+1,n-1,1)) * abs(e(i+1,n-1,1))^gamma1)/(DP(i+1)*h*q);
                        else
                            Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))- epsilon * ur - delta/q * (sign(S(i,n))+S(i,n) )) / (DP(i)*h)-dd2_hat(i,n)-dd11_hat(i,n);
                        end
                        
                        if flagCtrlSat
                            Ui = max(min(Ui, aM), am);
                            
                        end
                    end
                    
                end
            else%线性反馈控制
                
                ur = -Ki'*C*Xi_error_sum;
                if intOfUr(i) == 0
                    intOfUr(i) = Xi(2);
                else
                    intOfUr(i) = intOfUr(i) + ur * deltaT;
                end
                s(i,n) = Xi(2) - intOfUr(i);
                if flagCtrlSgn2Sat
                    Ui = ur;
                else
                    Ui = ur ;%去掉了这里面的滑模，需要加上的话粘贴到ur后面：- ks * sign(s(i,n))
                    
                end
            end
            
            if flagParaMismatch
                Ui = ctrlPara(1,i) * Ui + ctrlPara(2,i) * Xi(2)^2 + ctrlPara(3,i) * g;
            end
            if flagCtrlSat
                Ui = max(min(Ui, aM), am);
            end
        end
        
        x_wave_neighbor_average(i,n,:)=Xi_error_sum/DP(i);
        u_neighbor_average(i,n) = Ui_sum/DP(i);
        
        Gi=reshape(G(i,:,:),3,3);
        Fi=reshape(F(i,:,:),3,1);
        Ni=reshape(N(i,:,:),3,2);
        DDi=reshape(dd(i,n,:),2,1);
        NXi=Gi*Xi+Fi*Ui+Ni*DDi;
        if arvTimeCtrl(i) == 0 && NXi(1) >= -Dmz
            arvTimeCtrl(i) = n*Tstep;
        end
        
        u(i,n)=Ui;
        x(i,n,:)=NXi;
        %         x_error(i,n,:)=Xi_error;
    end
end

%% 绘图
p0=x0(1,:);
v0=x0(2,:);
a0=x0(3,:);
p=reshape(x(:,:,1),Nveh,Nstep);
v=reshape(x(:,:,2),Nveh,Nstep);
a=reshape(x(:,:,3),Nveh,Nstep);
dd1=reshape(dd(:,:,1),Nveh,Nstep);
dd2=reshape(dd(:,:,2),Nveh,Nstep);
v_r=v+dd1;

vv=[v0;v];
ddvv=vv(2:end,:)-vv(1:end-1,:);

t=(1:Nstep-tstart)*Tstep;
%t=t(tstart+1:end);
dp=p-kron(ones(Nveh,1),p0)+kron(ones(1,Nstep),kron((1:Nveh)',d0));
dv=v-kron(ones(Nveh,1),v0);
da=a-kron(ones(Nveh,1),a0);%位置速度加速度误差
DOdd1=dd1-dd1_hat;
DOdd2=dd2-dd2_hat;
DOdp=p-p_hat;
DOdv=v-v_hat;
% x3_error=reshape(x_error(3,:,:),Nstep,3);
space=reshape(e(:,:,1),Nveh,Nstep);
ddv=reshape(e(:,:,2),Nveh,Nstep);

if max(max(p))>2000
    doPic=true;
end

if doPic == true
    
    hf(1)=figure(1);
    hold on
    hp=plot(t,p(:,tstart+1:end),'linewidth',0.8);
    plot(t,p0(:,tstart+1:end),'k:','linewidth',0.8);
    legend('1','2','3','4','5','6','0','linewidth',0.3,'fontsize',4,'linewidth',0.3,'fontsize',4,'Location','southeast');  legend('boxoff');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{p_i }\rm[m]');
    xlim([0 max(t)])  ; ylim([0 1500]);
    set(gcf,'Position',[850 150 165 140]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%设置坐标轴的标签字号大小
    ylabel('\it{p_i }\rm[m]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
   % set(gca,'YTick',[15:2:25]) %? %改变x轴坐标间隔显示
    set(gca,'XTick',[0:10:55]) %? %改变x轴坐标间隔显示
  %  ,'FontWeight','bold'  字体加粗
    ha(1)=gca; 
%     hf(2)=figure(2);
%     hold on
%     hp=plot(t,dp(:,tstart+1:end));
%     %     legend(hp(1:2:end),'1','3','5','7','9');
%     legend(hp(1:1:end),'Veh 1','Veh 2','Veh 3','Veh 4','Veh 5','Veh 6');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{\Deltap}\rm[m]');
%     axes('Position',[0.18,0.62,0.28,0.25]); % 生成子图
%     plot(t,dp(:,tstart+1:end)); % 绘制局部曲线图
%     xlim([min(8),max(15)]); % 设置坐标轴范围
%     ha(2)=gca;
    
    hf(8)=figure(8);
    hold on
    hp=plot(t,DOdd1(:,tstart+1:end),'linewidth',0.8);
    legend('1','2','3','4','5','6','linewidth',0.3,'fontsize',4,'linewidth',0.3,'fontsize',4,'Location','southeast');  legend('boxoff');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{\omega_{i1}}\rm[m/s]');
    xlim([0 max(t)])  ; ylim([-0.2 0.2]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%设置坐标轴的标签字号大小
    ylabel('\it{\omega_{i1}}\rm[m/s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
   % set(gca,'YTick',[15:2:25]) %? %改变x轴坐标间隔显示
    set(gca,'XTick',[0:10:55]) %? %改变x轴坐标间隔显示
    axes('Position',[0.32,0.52,0.38,0.25]); % 生成子图
    plot(t,DOdd1(:,tstart+1:end)); % 绘制局部曲线图
    xlim([min(15),max(25)]); % 设置坐标轴范围
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    set(gcf,'Position',[850 150 165 140]);
  %  ,'FontWeight','bold'  字体加粗
    ha(8)=gca;
    
      hf(9)=figure(9);
    hold on
    hp=plot(t,DOdd2(:,tstart+1:end));
    legend('1','2','3','4','5','6','linewidth',0.3,'fontsize',4,'linewidth',0.3,'fontsize',4,'Location','southeast');  legend('boxoff');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{\omega_{i1}}\rm[m/s]');
    xlim([0 max(t)])  ; ylim([-0.5 0.5]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%设置坐标轴的标签字号大小
    ylabel('\it{\omega_{i2}}\rm[m/s^2]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
    set(gca,'YTick',[-0.5:0.2:0.5]) %? %改变x轴坐标间隔显示
    set(gca,'XTick',[0:10:55]) %? %改变x轴坐标间隔显示
    axes('Position',[0.32,0.52,0.38,0.25]); % 生成子图
    plot(t,DOdd2(:,tstart+1:end)); % 绘制局部曲线图
    xlim([min(15),max(25)]); % 设置坐标轴范围
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    set(gcf,'Position',[850 150 165 140]);
  %  ,'FontWeight','bold'  字体加粗
    ha(9)=gca;
    
%     hf(9)=figure(9);
%     hold on
%     hp=plot(t,DOdd2(:,tstart+1:end));
%     ylim([min(-0.5),max(0.5)]); % 设置坐标轴范围
%     legend(hp(1:1:end),'Veh 1','Veh 2','Veh 3','Veh 4','Veh 5','Veh 6');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{\omega_{i2}}\rm[m/s^2]');
%     axes('Position',[0.18,0.62,0.28,0.25]); % 生成子图
%     plot(t,DOdd2(:,tstart+1:end)); % 绘制局部曲线图
%     xlim([min(6),max(12)]); % 设置坐标轴范围
%     ha(9)=gca;
    
    %     hf(10)=figure(10);
    %     hold on
    %     hp=plot(t,s(:,tstart+1:end));
    %     legend(hp(1:1:end),'Veh 1','Veh 2','Veh 3','Veh 4','Veh 5','Veh 6');
    %     xlabel('\it{t}\rm[s]');
    %     ylabel('\it{s}\rm[m]');
    %     ha(10)=gca;
    
    hf(11)=figure(11);
    hold on
    hp=plot(t,u(:,tstart+1:end));
    legend(hp(1:1:end),'Veh 1','Veh 2','Veh 3','Veh 4','Veh 5','Veh 6');
    xlabel('\it{t}\rm[s]');
    ylabel('\it{u}\rm[m/s^2]');
    ha(11)=gca;
    
    hf(12)=figure(12);
    hold on
    hp=plot(t,v_r(:,tstart+1:end),'linewidth',0.8);
    plot(t,v0(:,tstart+1:end),'k:','linewidth',0.8);
    legend('0','linewidth',0.3,'fontsize',4);  legend('boxoff'); % 
    legend('1','2','3','4','5','6','0','linewidth',0.3,'fontsize',4);  legend('boxoff');
    xlabel('\it{t}\rm[s]');
    ylabel('\it{v_{ri} }\rm[m/s]');
    %     axes('Position',[0.32,0.52,0.38,0.25]); % 生成子图
    %     plot(t,v_r); % 绘制局部曲线图
    %     xlim([min(28),max(35)]); % 设置坐标轴范围
        xlim([0 max(t)])  ; ylim([15 25]);
    set(gcf,'Position',[850 150 165 140]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%设置坐标轴的标签字号大小
    ylabel('\it{v_{ri} }\rm[m/s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
    set(gca,'YTick',[15:2:25]) %? %改变x轴坐标间隔显示
  %  set(gca,'XTick',[0:3:15]) %? %改变x轴坐标间隔显示
  %  ,'FontWeight','bold'  字体加粗
    ha(12)=gca;
%     hf(13)=figure(13);
%     hold on
%     hp=plot(t,S(:,tstart+1:end));
%     legend(hp(1:1:end),'Veh 1','Veh 2','Veh 3','Veh 4','Veh 5','Veh 6');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{S}\rm[m]');
%     ha(13)=gca;
    
    hf(14)=figure(14);
    hold on
    hp=plot(t,space(:,tstart+1:end),'linewidth',0.7);
    %   ylim([min(-0.6),max(0.4)]); % 设置坐标轴范围
    %     legend(hp(1:2:end),'1','3','5','7','9');
    legend(hp(1:1:end),'1','2','3','4','5','6','linewidth',0.3,'fontsize',4);  legend('boxoff');
   % xlabel('\it{t}\rm[s]'); ylabel('\it{e_i }\rm[m]');
    xlim([0 max(t)])  ; ylim([-0.2 0.1]);  %ylim([-0.5 1]);
        set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%设置坐标轴的标签字号大小
    ylabel('\it{e_i }\rm[m]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
    axes('Position',[0.32,0.52,0.38,0.25]); % 生成子图
    plot(t,space(:,tstart+1:end)); % 绘制局部曲线图
    xlim([19 22])  ; ylim([-0.06 0.06]);
    set(gcf,'Position',[850 150 165 140]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on

  %  set(gca,'XTick',[0:3:15]) %? %改变x轴坐标间隔显示
    ha(14)=gca;
    
    hf(16)=figure(16);
    hold on
    hp=plot(t,ddvv(:,tstart+1:end),'linewidth',0.7);
    %  ylim([min(-0.8),max(0.5)]); % 设置坐标轴范围
    %     legend(hp(1:2:end),'1','3','5','7','9');
    legend(hp(1:1:end),'1','2','3','4','5','6','linewidth',0.3,'fontsize',4,'Location','southeast');legend('boxoff');
 %   xlabel('\it{t}\rm[s]');ylabel('\it{e_{v,i} }\rm[m/s]');
    xlim([0 max(t)])  ; ylim([-1.2 0.8]); %ylim([-1.5 1.2]);
    set(gcf,'Position',[850 150 165 140]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%设置坐标轴的标签字号大小
    ylabel('\it{e_{v,i} }\rm[m/s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
  %  set(gca,'XTick',[0:3:15]) %? %改变x轴坐标间隔显示
  %  ,'FontWeight','bold'  字体加粗
    ha(16)=gca;
%     figure;    参考绘图设置
% plot(t, Postion(:,1)-(x0(1:Time_sim/Tim_step,:)- d),'r','linewidth',1.5);hold on;
% plot(t, Postion(:,2)-(x0(1:Time_sim/Tim_step,:)- 2*d) ,'b','linewidth',1.5);hold on;
% plot(t, Postion(:,3)-(x0(1:Time_sim/Tim_step,:)- 3*d),'k','linewidth',1.5);hold on;
% plot(t, Postion(:,4)-(x0(1:Time_sim/Tim_step,:)- 4*d),'m','linewidth',1.5);hold on;
% plot(t, Postion(:,5)-(x0(1:Time_sim/Tim_step,:)- 5*d),'g','linewidth',1.5);hold on;
% % plot(t, Postion(:,6)-(x0(1:Time_sim/Tim_step,:)- 6*d),'r--','linewidth',2);hold on;
% % plot(t, Postion(:,7)-(x0(1:Time_sim/Tim_step,:)- 7*d),'b--','linewidth',2);hold on;
% h = legend('1','2','3','4','5','6','7');
% set(h,'box','off'); box off;
% xlabel('Time (s)');ylabel('Spacing error (m)');
% xlim([0 T])  ;ylim([-0.3 0.25])
% % axes('Position',[0.32,0.52,0.38,0.25]); % 生成子图
% % plot(t,DOdd1(:,tstart+1:end)); % 绘制局部曲线图
% % xlim([min(6),max(12)]); % 设置坐标轴范围
% set(gcf,'Position',[850 150 260 250]);
% set(gca,'linewidth',1,'fontsize',10,'fontname','times new roman');
%  % grid on
% xlabel('Time (s)','fontname', 'times new roman','fontSize',10);%设置坐标轴的标签字号大小
% ylabel('Spacing error (m)','fontname', 'times new roman','fontSize',10);
% set(gca,'XTick',[0:3:15]) %? %改变x轴坐标间隔显示 
end


%% Measurement of Effectiveness
maxep= max(max(abs(e(:,tstart:end,1))));
maxev= max(max(abs(e(:,tstart:end,2))));
aveep=(sum(sum(abs(e(:,tstart:end,1)))))/numel(e(:,tstart:end,1));
aveev=(sum(sum(abs(e(:,tstart:end,2)))))/numel(e(:,tstart:end,2));
[maxep,maxev,aveep,aveev]