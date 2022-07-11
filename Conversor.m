function [dil, dvc] = Conversor(d,il,x)
    %% Parametro
    Vin = 20; %Tensão de entrada
    Po = 10; %Potência da CPL
    R = 4; % Resistência de carga
    L = 1e-3; %Indutância
    rl = 0; % Resistência de enrolamento
    C = 2.2e-3; %Capacitância

    k1 = -1/L;
    k2 = rl*k1;
    k3 = -Vin*k1;
    k4 = 1/C;
    k5 = -k4/R;
    k6 = k4*Po;
    
    dil = k2*il + k1*x + d*k3;
    dvc = k4*il + k5*x + k6/x;
    
end