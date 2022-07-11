close all, clear all,clc
%% Parâmetros
Vin = 20; %Tensão de entrada
Vout = 12; %Tensão de saída
d = 0.6; % Ciclo de trabalho
PoV = [0, 10, 15, 30, 70]; %Potência da CPL
RV = [4,2,6,8,10 ]; % Resistência de carga
L = 10^-3; %Indutância
rl = 0; % Resistência de enrolamento
CV = [2.2,0.5, 1.0, 1.5,  2.5]; %Capacitância
vcin = 0.01;
ilin = 0.01;
u0 = 0.05;
v=0;
flag = 0;
for j = 1:1
    about = 'Variando ';
    for k=1:1
        R =4; Po = 10; C = 2.2*10^-3;
        if flag == 0
            R = RV(k);
            about = ['Variando R = ' int2str(R) ' Ohms'];
        elseif flag == 1
            Po = PoV(k);
            about = ['Variando Po = ' int2str(Po) ' W'];
        elseif flag == 2
            C = CV(k)*10^-3;
            about = ['Variando C = ' num2str(C) ' F'];
        end
        vc0 = d*Vin; % -rl*il = 0 já q rl = 0;
        il0 = vc0/R - Po/vc0;
        

        k1 = -1/L;
        k2 = rl*k1;
        k3 = -Vin*k1;
        k4 = 1/C;
        k5 = -k4/R;
        k6 = k4*Po;

        %% Discretização por Euler
        Ts = 1e-5; % Tempo de amostragem
        T = 0.2; % Tempo de simulação
        t = 0:Ts:T; % Vetor do tempo
        P0 = (T/2)/Ts+1; % Step time
        di = zeros(1,length(t)); % Entrada
        di(1:end)=d;
        di(P0:end) = d+u0;

        vc = vcin*ones(1,length(t));
        vc_1 = vcin; % Altura passada
        vc_0 = 0; % Altura atual
        il = ilin*ones(1,length(t));
        il_1 = ilin; % Altura passada
        il_0 = 0; % Altura atual

        d_1 = d; % Entrada passada
        d_0 = d; % Entrada autal


        for ii=1:length(t)-1
            il_0 = il_1 + Ts*(k2*il_1 + k1*vc_1+k3*d_1);     
            vc_0 = vc_1 + Ts*k4*il_1+Ts*k5*vc_1 + Ts*k6*(1/vc_1);
            il_1 = il_0;
            vc_1 = vc_0;
            vc(ii+1) = vc_0;
            d_1 = d_0;
            d_0 = di(ii+1);
        end

        %% Método de Rugen-Kutta de 4ª ordem
        x = vcin*ones(1,length(t));
        z = ilin*ones(1,length(t));
        u = di;
        for ii=1:length(t)-1
            u_1 = u(ii);
            [i_1, k_1] = Conversor(u_1,z(ii),x(ii));
            u_2 = u(ii);
            [i_2, k_2] = Conversor(u_2,z(ii)+Ts*i_1/2,x(ii)+Ts*k_1/2);
            u_3 = u(ii);
            [i_3, k_3] = Conversor(u_3,z(ii)+Ts*i_2/2,x(ii)+Ts*k_2/2);
            u_4 = u(ii);
            [i_4, k_4] = Conversor(u_4,z(ii)+Ts*i_3,x(ii)+Ts*k_3);
            x(ii+1) = x(ii) + Ts*(k_1+2*k_2+2*k_3+k_4)/6;
            z(ii+1) = z(ii) + Ts*(i_1+2*i_2+2*i_3+i_4)/6;
        end
        figure("Name", ['Modelo não-linearizado e linearizado ' about])
        outsim = sim('Lab1Exp1Passo2',T);
        subplot(3,1,1);
        plot(t, vc, t, x, 'g', t, outsim.DiagramaBlocos(:,2));
        legend('Discretização por Euler', 'Método de Rugen-Kutta de 4ª ordem', 'Diagrama de Blocos');

        A = [k2 k1;k4 k5-k4*1/vc0^2];
        B = [Vin/L d/L 0; 0 0 1/(C*vc0)];
        B1 = [Vin/L; 0];
        C1 = [0 1];
        D = [0];
        x0 = [ilin vcin];
        sys = ss(A,B1,C1,D);
        sysTf = tf(sys);
        u2 = Vin*ones(1,length(t));
        u3 = Po*ones(1,length(t));
        ut = [u; u2; u3];
        [yss,t,xss] = lsim(sys,u,t,x0);
        subplot(3,1,2);
        plot(t, yss);
        xlim([0 T]);
        [yft,t,xft] = lsim(sysTf,u,t);
        subplot(3,1,3);
        plot(t, yft);
        xlim([0 T]);
        
        
    end
    flag = flag+1;
end