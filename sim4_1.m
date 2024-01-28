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
% ���ڵİ汾����DO�������Ҳ����������ȶ�                 qin���������������Ҫ��������Ӧ����
%  �Լ��Ŀ�����������DO�Աȣ�ͬʱ�������������Աȣ������������Ŷ���Ч�����У�   DO���������ٵ�һ��
%% ���̿���
doPic = true;
doFeedforward = false; %�Ƿ�ǰ��
isDelayed = false; %ǰ���е�������������Ƿ���ʱ��
isRandom = 0 ; %1�������ʼ״̬,0����ƫ���ʼ״̬
isPseudoRandom = false; %1���̶��������ʼ״̬��0���ǹ̶��������ʼ״̬
isHeterogeneous = true; %����ָ���Ƿ�����
doDesiredAccelerationInput = false; %�캽�����з�ʽ��trueΪ�������ٶ����룬falseΪʵ�ʼ��ٶ�����
doAccelerate = false; %�캽���Ƿ�ʼ�ռ���
IndexTopo = 1 ; % �������� 1.PF 2.PLF 3.TPF 4.TPLF 5.MF
flagParaMismatch = 0; % �Ƿ���ڲ�����ƥ��
flagCtrlSgn2Sat = 0  ; % sat ��� sgn
flagCtrlSat = 1 ; % �Ƿ�������뱥��
flagmisdis = 1 ;  %�Ƿ����ʧ���Ŷ�dd1,dd2
flagDO = 1  ;   %�Ƿ�����Ŷ��۲���
flagfiniteDO = 1 ; %   1����ʱ��DO,  0����wangDO��ͬʱ�л�Ϊwang2020������
flagFinitefTimeCtrl = 1 ; % ��������ʱ����ƣ�0Ϊ���Է������ƣ�qin³����������
%% ��ʼ��
ddes= 8 ;% ��С��ȫ���롣������6
v0= 16 ; %��ʼ����   6
v00 = 16;%�캽����ʼ���� 6
D0=[ddes;0;0];%״̬ƫ��
h= 0.5 ;%%ʱ��
d0=ddes+h*v0; %��ʼ����  12��������9
C=diag([1,1,1]);
k11 = 2.5; k22 = 1.26; k33 = 0.23 ;   %%%   ���Է�����������
% k11 = 0.1581; k22 = 0.5142; k33 = 0.2571 ;   %%%   ���Է�����������  ���¼���
% k11 = 0.84; k22 = 2.41; k33 = 1.16 ;   %%%   ���Է�����������  commen
% k11 = 1; k22 = 0.5; k33 = 0;
k1 = 0.15 ;   % ��Ӧ�����е�k_i   0.6
epsilon = 1   ; %   0.3  0.5  �������������û�У�����������������ȡΪ1
delta = 0.5 ; %    ��Ӧ�����е�gamma  0.45
q= 0.8 ;%��ϻ�ģ����0.7����Ӧ�����е�betta
ppp = 0.3 ;  %% 0.4  ��Ӧ�����е�alfa_iS
% D0 = 10;
g = 9.81;
Dmz = 5;
Dcz = 150;
am = -3; %minimum acceleration
aM = 3; %maximum acceleration
lanmda10=  10  ;  %%ԭ 10   5
lanmda11=  8  ;%%%ԭ 8     6
lanmda12=  7  ;%% ԭ 7     4
lanmda20= 5 ; %ԭ  5     6    (((5   5
lanmda21= 6 ; % ԭ 6     8
lip1= 0.5 ;  %0.35
lip2= 1.5 ;  %% 0.6
%ppp = 0.3;   %%ԭ gamma2*gamma3/(2*gamma3 - gamma2)  wang2020��������Ӧ����ppp
tstart= 1000; %��ͼ��ʼʱ��
tdd1 = 20+tstart/100  ;   %�Ŷ�1����ʱ��
tdd2= 20+tstart/100   ;
fdd1= 0.3 ;
fdd2= 1.5 ;  %%%�Ŷ���ֵ
tdd3 = 20+tstart/100  ;   %�Ŷ�3����ʱ��
tdd4= 20+tstart/100   ;
fdd3= 0 ;
fdd4= 0 ;  %%%�Ŷ���ֵ
l=0.1;  %wang DO����

Nveh= 6 ; %������
if Nveh<3
    msgbox('Nveh>=3')
end

Tstep=0.01; %���沽��
deltaT = Tstep;
Nstep=tstart+5500; %���沽��
tau = 0.3;


K=zeros(3,Nveh);%��������
G=zeros(Nveh,3,3);%��ɢ״̬����
F=zeros(Nveh,3,1);%��ɢ�������
x=zeros(Nveh,Nstep,3);%����״̬
N=zeros(Nveh,3,2);%�Ŷ�ϵ��
dd=zeros(Nveh,Nstep,2);%��ɢ�Ŷ�
s=zeros(Nveh,Nstep);%��ɢ��ģ
S=zeros(Nveh,Nstep);%��ɢ��ϻ�ģ
x_wave_neighbor_average=zeros(Nveh,Nstep,3);%����״̬���
u_neighbor_average=zeros(Nveh,Nstep); %��������
e=zeros(Nveh,Nstep,2);%�������
Vi_sum=zeros(Nveh,Nstep); %�ڳ��ٶ�֮��
%����Ϊ�۲������ֲ���
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


z=zeros(Nveh,Nstep);% wang DO �м���  Ҳ���������е�z_hat
z_dot=zeros(Nveh,Nstep);

for i=2:Nveh
    p_hat(i,1)=p_hat(i-1,1)-d0;
end
v_hat(:,1)=v_hat(:,1)+v0;%�۲�����ʼ��

%% ����
M=zeros(Nveh,Nveh); %�ڽӾ���
P=zeros(Nveh,Nveh); %ǣ������

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

MP=M+P; %ͼ����������
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

randp = [0.3192    0.3129   -0.8649   -0.0301   -0.1649    0.6277    1.0933]; % α���λ���Ŷ�
randv = [1.4384    0.3252   -0.7549    1.3703   -1.7115   -0.1022   -0.2414]; % α����ٶ��Ŷ�
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
            x(i,1,:)=[randp(i)-d0*i; v0+randv(i) ;0];%��ʼ״̬
        else
            x(i,1,:)=[normrnd(0,1)-d0*i; v0+normrnd(0,1) ;0];%��ʼ״̬
        end
    else
        x(i,1,:)=[-d0*i; v0;0];%��ʼ״̬
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
x0=zeros(3,Nstep); %�캽��״̬
x0(:,1)=[d0*Nveh+h*v0;v00;0];%�캽��״̬��ʼ������һ�У�
u0=zeros(1,Nstep);
u=zeros(Nveh,Nstep);

%% �ջ�����ѧ
for n=2:Nstep %����
    x0(:,n)=x0(:,n-1);%�캽��״̬����
    
    %����Ϊ�۲�������
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
    if  n<=4/Tstep+tstart   %  ԭ  1
        x0(3,n)=0;
        u0(n)=0;
    else
        if doAccelerate %�㶨����
            x0(3,n)=1;
            u0(n)=1;
        else %���١�����
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
    if doDesiredAccelerationInput %�캽�������������ٶ�����
        x0(:,n)=Gi*x0(:,n-1)+Fi*u0(n);%�캽��״̬����
    else %�캽������ʵ�ʼ��ٶ�����
        x0(2,n)=x0(2,n)+x0(3,n)*Tstep;%�캽���ٶȸ���
        x0(1,n)=x0(1,n)+x0(2,n)*Tstep+1/2*x0(3,n)*Tstep^2;%�캽��λ�ø���
    end
    if arvTimeOptm(1) == 0 && x0(1,n) >= -Dmz
        arvTimeOptm(1) = n*Tstep;
        arvTimeOptm = arvTimeOptm(1) + (1:Nveh)'*d0/v00;
    end
    
    for i=1:Nveh %MP�������,����n����
        Xi=reshape(x(i,n-1,:),3,1);%�Գ�״̬
        Xi_neighbor=zeros(3,1); %�ڳ�״̬
        Xi_error_sum=zeros(3,1); %״̬���
        Ui_neighbor=0; %�ڳ�����
        Ui_sum=0; %�ڳ�����֮��
        Vi_neighbor=0; %�ڳ��ٶ�
        
        if flagmisdis
            dd(i,n,1)=  dd(i,n,1)+fdd1*sin(1*n*Tstep)*exp(-(n*Tstep-tdd1)^2) +fdd3*sin(0.5*n*Tstep) ;%%  dd(i,n,1)+0.25*sin(n*Tstep)
            dd(i,n,2)= dd(i,n,2)+fdd2*sin(3*n*Tstep)*exp(-(n*Tstep-tdd2)^2)+fdd4*sin(0.5*n*Tstep) ;%�Ŷ�����   dd(i,n,2)+0.2*sin(n*Tstep)   -0.2*i
        else
            dd(i,n,1)=  0  ;  %%  dd(i,n,1)+0.25*sin(n*Tstep)
            dd(i,n,2)= 0 ;
        end
        
        for j=1:i %MP�������,����n���ڳ�  ԭ j=1:Nveh
            if MP(i,j)~=0 %�����ڳ�
                if i==j %�ڳ�Ϊ�캽��
                    Xi_neighbor=x0(:,n);
                    Xi_neighbor(1)=Xi_neighbor(1)-h*(x(i,n-1,2)+dd(i,n,1))-Tstep*v0;%%  %  -Tstep*v0��Ϊ��������ɢ�����µĳ�ʼ���
                    Xi_neighbor(2)=Xi_neighbor(2);%��������Ѹ���
                    Vi_neighbor=x0(2,n);
                    if doDesiredAccelerationInput %�캽�������������ٶ�����
                        Ui_neighbor=u0(n);
                        Ui_neighbor_delayed=u0(n-1);
                    else %�캽������ʵ�ʼ��ٶ�����
                        Ui_neighbor=x0(3,n);
                        Ui_neighbor_delayed=x0(3,n-1);
                    end
                else %�ڳ����캽��
                    Xi_neighbor=reshape(x(j,n-1,:),3,1);
                    %                   if i==1
                    %                      Xi_neighbor(1)=0;
                    %                   else
                    Xi_neighbor(1)=Xi_neighbor(1)-h*(x(i,n-1,2)+dd(i,n,1));
                    %                   end
                    Xi_neighbor(2)=Xi_neighbor(2) ;%��������Ѹ���
                    Ui_neighbor=u(j,n)+dd1_hat_dot(j,n);%�Ѹ��ݿ���������
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
            else %�������ڳ�
                continue;
            end
        end
        Ki=K(:,i);
        e(i,n,1) = Xi_error_sum(1);
        e(i,n,2) = Xi_error_sum(2);
        
        %����۲���
        if flagDO
            if  flagfiniteDO   %����do
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
            else  %����do
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
            if flagFinitefTimeCtrl%����ʱ�����
                
                sumOfNghbAcc = Ui_sum;
                ur =  k1 * e(i,n,1)  ;   %�޸�   ���������ʣ�����
                %                 if intOfUr(i) == 0
                %                     intOfUr(i) = 0 ;%%e(i,n,2)/epsilon
                %                 else
                intOfUr(i) = intOfUr(i) + ur * deltaT;%ur�Ļ���
                %                 end
                s(i,n) = e(i,n,1) + epsilon * intOfUr(i);
                if i<Nveh
                    S(i,n) = q*s(i,n)-s(i+1,n-1);
                else
                    S(i,n) = q*s(i,n);
                end
                if  flagfiniteDO   %����do  �����ﲻֹ�л�DO�����л���������wang2020��
                    if flagCtrlSgn2Sat
                        Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))+ epsilon * ur - delta/q * sign(S(i,n)) * abs(S(i,n))^ppp ) / (DP(i)*h)-dd2_hat(i,n);
                    else
                        if i<Nveh
                            if n<30        Vi_sum(i+1,n-1)=x(i+1,n-1,2)+dd1_hat(i+1,n); else  ;  end
                            Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))- epsilon * ur - delta/q * sign(S(i,n)) * abs(S(i,n))^ppp ) / (DP(i)*h)-dd2_hat(i,n)-dd11_hat(i,n) ...
                                +(DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))+k1 * e(i+1,n-1,1))/(DP(i+1)*h*q);% ȥ����e_i+1�ĵ�����һ��  DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))+
                            %+(DP(i+1)*(x(i+1,n-1,2)+dd1_hat(i+1,n))-Vi_sum(i+1,n-1)+DP(i+1)*h*(u(i+1,n-1)+dd2_hat(i+1,n))- k1 * sign(e(i+1,n-1,1)) * abs(e(i+1,n-1,1))^gamma1)/(DP(i+1)*h*q);
                        else
                            Ui = ( Vi_sum(i,n)-DP(i)*(x(i,n-1,2)+dd1_hat(i,n))- epsilon * ur - delta/q * sign(S(i,n)) * abs(S(i,n))^ppp ) / (DP(i)*h)-dd2_hat(i,n)-dd11_hat(i,n);
                        end
                        
                        if flagCtrlSat
                            Ui = max(min(Ui, aM), am);
                        end
                    end
                else   %%  ����do  wang2019������
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
            else%���Է�������
                
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
                    Ui = ur ;%ȥ����������Ļ�ģ����Ҫ���ϵĻ�ճ����ur���棺- ks * sign(s(i,n))
                    
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

%% ��ͼ
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
da=a-kron(ones(Nveh,1),a0);%λ���ٶȼ��ٶ����
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
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%����������ı�ǩ�ֺŴ�С
    ylabel('\it{p_i }\rm[m]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
   % set(gca,'YTick',[15:2:25]) %? %�ı�x����������ʾ
    set(gca,'XTick',[0:10:55]) %? %�ı�x����������ʾ
  %  ,'FontWeight','bold'  ����Ӵ�
    ha(1)=gca; 
%     hf(2)=figure(2);
%     hold on
%     hp=plot(t,dp(:,tstart+1:end));
%     %     legend(hp(1:2:end),'1','3','5','7','9');
%     legend(hp(1:1:end),'Veh 1','Veh 2','Veh 3','Veh 4','Veh 5','Veh 6');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{\Deltap}\rm[m]');
%     axes('Position',[0.18,0.62,0.28,0.25]); % ������ͼ
%     plot(t,dp(:,tstart+1:end)); % ���ƾֲ�����ͼ
%     xlim([min(8),max(15)]); % ���������᷶Χ
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
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%����������ı�ǩ�ֺŴ�С
    ylabel('\it{\omega_{i1}}\rm[m/s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
   % set(gca,'YTick',[15:2:25]) %? %�ı�x����������ʾ
    set(gca,'XTick',[0:10:55]) %? %�ı�x����������ʾ
    axes('Position',[0.32,0.52,0.38,0.25]); % ������ͼ
    plot(t,DOdd1(:,tstart+1:end)); % ���ƾֲ�����ͼ
    xlim([min(15),max(25)]); % ���������᷶Χ
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    set(gcf,'Position',[850 150 165 140]);
  %  ,'FontWeight','bold'  ����Ӵ�
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
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%����������ı�ǩ�ֺŴ�С
    ylabel('\it{\omega_{i2}}\rm[m/s^2]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
    set(gca,'YTick',[-0.5:0.2:0.5]) %? %�ı�x����������ʾ
    set(gca,'XTick',[0:10:55]) %? %�ı�x����������ʾ
    axes('Position',[0.32,0.52,0.38,0.25]); % ������ͼ
    plot(t,DOdd2(:,tstart+1:end)); % ���ƾֲ�����ͼ
    xlim([min(15),max(25)]); % ���������᷶Χ
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    set(gcf,'Position',[850 150 165 140]);
  %  ,'FontWeight','bold'  ����Ӵ�
    ha(9)=gca;
    
%     hf(9)=figure(9);
%     hold on
%     hp=plot(t,DOdd2(:,tstart+1:end));
%     ylim([min(-0.5),max(0.5)]); % ���������᷶Χ
%     legend(hp(1:1:end),'Veh 1','Veh 2','Veh 3','Veh 4','Veh 5','Veh 6');
%     xlabel('\it{t}\rm[s]');
%     ylabel('\it{\omega_{i2}}\rm[m/s^2]');
%     axes('Position',[0.18,0.62,0.28,0.25]); % ������ͼ
%     plot(t,DOdd2(:,tstart+1:end)); % ���ƾֲ�����ͼ
%     xlim([min(6),max(12)]); % ���������᷶Χ
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
    %     axes('Position',[0.32,0.52,0.38,0.25]); % ������ͼ
    %     plot(t,v_r); % ���ƾֲ�����ͼ
    %     xlim([min(28),max(35)]); % ���������᷶Χ
        xlim([0 max(t)])  ; ylim([15 25]);
    set(gcf,'Position',[850 150 165 140]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%����������ı�ǩ�ֺŴ�С
    ylabel('\it{v_{ri} }\rm[m/s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
    set(gca,'YTick',[15:2:25]) %? %�ı�x����������ʾ
  %  set(gca,'XTick',[0:3:15]) %? %�ı�x����������ʾ
  %  ,'FontWeight','bold'  ����Ӵ�
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
    %   ylim([min(-0.6),max(0.4)]); % ���������᷶Χ
    %     legend(hp(1:2:end),'1','3','5','7','9');
    legend(hp(1:1:end),'1','2','3','4','5','6','linewidth',0.3,'fontsize',4);  legend('boxoff');
   % xlabel('\it{t}\rm[s]'); ylabel('\it{e_i }\rm[m]');
    xlim([0 max(t)])  ; ylim([-0.2 0.1]);  %ylim([-0.5 1]);
        set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%����������ı�ǩ�ֺŴ�С
    ylabel('\it{e_i }\rm[m]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
    axes('Position',[0.32,0.52,0.38,0.25]); % ������ͼ
    plot(t,space(:,tstart+1:end)); % ���ƾֲ�����ͼ
    xlim([19 22])  ; ylim([-0.06 0.06]);
    set(gcf,'Position',[850 150 165 140]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on

  %  set(gca,'XTick',[0:3:15]) %? %�ı�x����������ʾ
    ha(14)=gca;
    
    hf(16)=figure(16);
    hold on
    hp=plot(t,ddvv(:,tstart+1:end),'linewidth',0.7);
    %  ylim([min(-0.8),max(0.5)]); % ���������᷶Χ
    %     legend(hp(1:2:end),'1','3','5','7','9');
    legend(hp(1:1:end),'1','2','3','4','5','6','linewidth',0.3,'fontsize',4,'Location','southeast');legend('boxoff');
 %   xlabel('\it{t}\rm[s]');ylabel('\it{e_{v,i} }\rm[m/s]');
    xlim([0 max(t)])  ; ylim([-1.2 0.8]); %ylim([-1.5 1.2]);
    set(gcf,'Position',[850 150 165 140]);
    set(gca,'linewidth',0.5,'fontsize',5,'fontname','times new roman','FontWeight','bold');
    % grid on
    xlabel('\it{t}\rm[s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');%����������ı�ǩ�ֺŴ�С
    ylabel('\it{e_{v,i} }\rm[m/s]','fontname', 'times new roman','fontSize',7,'FontWeight','bold');
  %  set(gca,'XTick',[0:3:15]) %? %�ı�x����������ʾ
  %  ,'FontWeight','bold'  ����Ӵ�
    ha(16)=gca;
%     figure;    �ο���ͼ����
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
% % axes('Position',[0.32,0.52,0.38,0.25]); % ������ͼ
% % plot(t,DOdd1(:,tstart+1:end)); % ���ƾֲ�����ͼ
% % xlim([min(6),max(12)]); % ���������᷶Χ
% set(gcf,'Position',[850 150 260 250]);
% set(gca,'linewidth',1,'fontsize',10,'fontname','times new roman');
%  % grid on
% xlabel('Time (s)','fontname', 'times new roman','fontSize',10);%����������ı�ǩ�ֺŴ�С
% ylabel('Spacing error (m)','fontname', 'times new roman','fontSize',10);
% set(gca,'XTick',[0:3:15]) %? %�ı�x����������ʾ 
end


%% Measurement of Effectiveness
maxep= max(max(abs(e(:,tstart:end,1))));
maxev= max(max(abs(e(:,tstart:end,2))));
aveep=(sum(sum(abs(e(:,tstart:end,1)))))/numel(e(:,tstart:end,1));
aveev=(sum(sum(abs(e(:,tstart:end,2)))))/numel(e(:,tstart:end,2));
[maxep,maxev,aveep,aveev]