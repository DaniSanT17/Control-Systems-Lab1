clc, close all, close
%% Parâmetros
Vin = 20; % V Tensão de entrada
Vout = 12; % V Tensão de saída
d = 0.6; % Ciclo de trabalho
Po = 10; % W Potência da CPL
R = 4; % Ohms Resistência de carga
L = 1e-3; % H Indutância
rl = 0; % Ohms Resistência de enrolamento
C = 2.2e-3; % F Capacitância
Vci = 1e-1;
ILi = 0;
u0 = d*0;

Ts = 1e-5; % Tempo de amostragem
T = 0.14; % Tempo de simulação
To = 0.08;
t = 0:Ts:T; % Vetor do tempo

il = 0;
Vc0 = -rl*il + d*Vin;
il0 = Vc0/R - Po/Vc0;

K1 = 1/L;
K2 = -K1*rl;
K3 = 1/C;
K4 = -K3/R;

out_sim = sim('Lab1Test', T);

%% Discretização por Euler
P0 = To/Ts+1; % Step time
di = zeros(1,length(t)); % Entrada
di(1:P0+1)=d; 
di(P0:end) = d+u0; % Acréscimo de 5% na entrada

vcE = Vci*ones(1,length(t)); % Vetor saída Vc
vc_1 = Vci; % Tensão passada
vc_0 = 0; % Tensão atual
il = ILi*ones(1,length(t)); % Vetor saída iL
il_1 = ILi; % Corrente passada
il_0 = 0; % Corrente atual

d_1 = d; % Entrada passada
d_0 = d; % Entrada autal


for ii=1:length(t)-1
    il_0 = il_1 + Ts*(K2*il_1 - K1*vc_1 + K1*Vin*d_1);     
    vc_0 = vc_1 + Ts*(K3*il_1 + K4*vc_1 + K3*Po/vc_1);
    il_1 = il_0;
    vc_1 = vc_0;
    vcE(ii+1) = vc_0;
    d_1 = d_0;
    d_0 = di(ii+1);
end

%% Método de Rugen-Kutta de 4ª ordem
vcRK = Vci*ones(1,length(t)); 
ilRK = ILi*ones(1,length(t));

for ii=1:length(t)-1
    u_1 = di(ii);
    [l_1, k_1] = CC(u_1,ilRK(ii),vcRK(ii));
    [l_2, k_2] = CC(u_1,ilRK(ii)+Ts*l_1/2,vcRK(ii)+Ts*k_1/2);
    [l_3, k_3] = CC(u_1,ilRK(ii)+Ts*l_2/2,vcRK(ii)+Ts*k_2/2);
    [l_4, k_4] = CC(u_1,ilRK(ii)+Ts*l_3,vcRK(ii)+Ts*k_3);
    
    ilRK(ii+1) = ilRK(ii) + Ts*(l_1+2*l_2+2*l_3+l_4)/6;
    vcRK(ii+1) = vcRK(ii) + Ts*(k_1+2*k_2+2*k_3+k_4)/6;
    
end


%% Espaço de Estados e Função de Transferência
Ass = [K2 -K1; K3 K4-K3*Po/(Vc0^2)];
Bss = [K1*Vin;  0];
Css = [0 1];
Dss = [0];

sys = ss(Ass,Bss,Css,Dss);
sysTf = tf(sys);
[yss,t,xss] = lsim(sys,di,t,[ILi Vci]);
[yft,t,xft] = lsim(sysTf,di,t);
figure();
plot(t,out_sim.Data(:,2),t,vcE,t, vcRK);
legend({'Diagrama de Blocos','Discretização por Euler', 'Rugen-Kutta de 4ª ordem'},'FontSize',14 );
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);

figure();
plot(t, yss, t, yft);
legend({'Espaço de Estados','Função de Transferência'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);


figure();
plot(t, yft,t, out_sim.Data(:,2));
legend({'Função de Transferência', 'Diagrama de Blocos'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);

%% Variando o R
Rvcl = [yft, yft, yft, yft];
Rvcnl = [out_sim.Data(:,2), out_sim.Data(:,2), out_sim.Data(:,2), out_sim.Data(:,2)];
Rs = [2,6,8,10];
for ii=1:4
    R = Rs(ii);
    Vc0 = d*Vin;
    il0 = Vc0/R - Po/Vc0;
    K1 = 1/L;
    K2 = -K1*rl;
    K3 = 1/C;
    K4 = -K3/R;
    out_sim = sim('Lab1Test', T);
    Rvcnl(:,ii) = out_sim.Data(:,2);
    
    Ass = [K2 -K1; K3 K4-K3*Po/(Vc0^2)];
    Bss = [K1*Vin;  0];
    Css = [0 1];
    Dss = [0];
    sys = ss(Ass,Bss,Css,Dss);
    sysTf = tf(sys)
    [yft,t,xft] = lsim(sysTf,di,t);
    Rvcl(:,ii) = yft;
    
end
figure();
plot(t, Rvcnl(:,1),t, Rvcnl(:,2), t, Rvcnl(:,3), t, Rvcnl(:,4));
legend({'Diagrama de Blocos R = 2 Ohms', 'Diagrama de Blocos R = 6 Ohms', 'Diagrama de Blocos R = 8 Ohms', 'Diagrama de Blocos R = 10 Ohms'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);
figure();
plot(t, Rvcl(:,1),t, Rvcl(:,2), t, Rvcl(:,3), t, Rvcl(:,4));
legend({'Função de Transferência R = 2 Ohms', 'Função de Transferência R = 6 Ohms', 'Função de Transferência R = 8 Ohms', 'Função de Transferência R = 10 Ohms'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);

%% Variando o Po
Rvcl = [yft, yft, yft, yft];
Rvcnl = [out_sim.Data(:,2), out_sim.Data(:,2), out_sim.Data(:,2), out_sim.Data(:,2)];
Ps = [0,15,30,70];
for ii=1:4
    R = 4;
    Po = Ps(ii);
    Vc0 = d*Vin;
    il0 = Vc0/R - Po/Vc0;
    K1 = 1/L;
    K2 = -K1*rl;
    K3 = 1/C;
    K4 = -K3/R;
    out_sim = sim('Lab1Test', T);
    Rvcnl(:,ii) = out_sim.Data(:,2);
    
    Ass = [K2 -K1; K3 K4-K3*Po/(Vc0^2)];
    Bss = [K1*Vin;  0];
    Css = [0 1];
    Dss = [0];
    sys = ss(Ass,Bss,Css,Dss);
    sysTf = tf(sys);
    [yft,t,xft] = lsim(sysTf,di,t);
    Rvcl(:,ii) = yft;
    
end
figure();
plot(t, Rvcnl(:,1),t, Rvcnl(:,2), t, Rvcnl(:,3), t, Rvcnl(:,4));
legend({'Diagrama de Blocos P = 0 W', 'Diagrama de Blocos P = 15 W', 'Diagrama de Blocos P = 30 W', 'Diagrama de Blocos P = 70 W'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);
figure();
plot(t, Rvcl(:,1),t, Rvcl(:,2), t, Rvcl(:,3), t, Rvcl(:,4));
legend({'Função de Transferência P = 0 W', 'Função de Transferência P = 15 W', 'Função de Transferência P = 30 W', 'Função de Transferência P = 70 W'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);

%% Variando o C
Rvcl = [yft, yft, yft, yft];
Rvcnl = [out_sim.Data(:,2), out_sim.Data(:,2), out_sim.Data(:,2), out_sim.Data(:,2)];
Cs = [0.5,1,1.5,2.5];
for ii=1:4
    Po = 10;
    C = Cs(ii)*10^-3;
    Vc0 = d*Vin;
    il0 = Vc0/R - Po/Vc0;
    K1 = 1/L;
    K2 = -K1*rl;
    K3 = 1/C;
    K4 = -K3/R;
    out_sim = sim('Lab1Test', T);
    Rvcnl(:,ii) = out_sim.Data(:,2);
    
    Ass = [K2 -K1; K3 K4-K3*Po/(Vc0^2)];
    Bss = [K1*Vin;  0];
    Css = [0 1];
    Dss = [0];
    sys = ss(Ass,Bss,Css,Dss);
    sysTf = tf(sys);
    [yft,t,xft] = lsim(sysTf,di,t);
    Rvcl(:,ii) = yft;
    
end
figure();
plot(t, Rvcnl(:,1),t, Rvcnl(:,2), t, Rvcnl(:,3), t, Rvcnl(:,4));
legend({'Diagrama de Blocos C = 0.5 mF', 'Diagrama de Blocos C = 1.0 mF', 'Diagrama de Blocos C = 1.5 mF', 'Diagrama de Blocos C = 2.5 mF'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);
figure();
plot(t, Rvcl(:,1),t, Rvcl(:,2), t, Rvcl(:,3), t, Rvcl(:,4));
legend({'Função de Transferência C = 0.5 mF', 'Função de Transferência C = 1.0 mF', 'Função de Transferência C = 1.5 mF', 'Função de Transferência C = 2.5 mF'}, 'FontSize',14);
xlabel('Tempo (s)', 'FontSize', 14);
ylabel('Tensão (V)', 'FontSize', 14);


%% Exp2
R = 4;
Po = 10;
C = 2.2e-3;

numFT = Vin/(C*L);
denC1 = 1;
denC2 = (1/(R*C)+Po/(C*Vc0^2)+rl/L);
denC3 = ((Po*rl/Vc0^2+rl/R+1)*1/(L*C));
H = 1;
ki = 1;
kp = 1;
kd = 1;
