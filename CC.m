function [dil, dvc] = CC(d,il,x)
    %% Parametro
    Vin = 20; %Tensão de entrada
    Po = 10; %Potência da CPL
    R = 4; % Resistência de carga
    L = 1e-3; %Indutância
    rl = 0; % Resistência de enrolamento
    C = 2.2e-3; %Capacitância

    K1 = 1/L;
    K2 = -K1*rl;
    K3 = 1/C;
    K4 = -K3/R;
    
    dil = K2*il - K1*x + d*K1*Vin;
    dvc = K3*il + K4*x + K3*Po/x;
    
end