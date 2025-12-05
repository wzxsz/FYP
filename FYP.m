clear
%% constant
Voc = 43.5;% Open circuit voltage V
Isc = 4.75;% Short circuit current A
Vmp = 34.5;% Maximum power voltage V
Imp = 4.35;% Maximum power current A
Ki = 0.00065 * Isc;% Temperature coefficient of ISC A/°C
Kv = -0.16;        % Temperature coefficient of VOC V/°C
q  = 1.602e-19;    % electron charge (C)
Kb = 1.380649e-23; % Boltzmann constant (J/K)
Ns = 72;% Number of series connected cells
Gr = 1000;  % Nominal irradiance level (W/m2) 
G = 1000;   % Operating irradiance level for first PV cell(W/m2)
G1 = 500;   % Operating irradiance level for second PV cell(W/m2)
Tr = 298.15;% Reference cell temperature(K) of 25°C
Tk = 298.15;% Operating temperature(K) suppose as 25°C 
Tdif = Tk - Tr;
A = 1.2;   % Diode ideality factor
Ego = 1.15;% Band gap energy for the silicon semiconductor eV
Rs = 1;    % Series connected resistance
Rp = 1e3; % Parallel connected resistance 
Vf  = 0.6; % Forward voltage for diode
Ron = 1e-2;% Forward resistance for diode
%% function
Iph_light = (Isc+Ki*Tdif) * G / Gr;
Iph_shaded = (Isc+Ki*Tdif) * G1 / Gr;
Irs = Isc / (exp(q*Voc/(Kb*A*Tk*Ns) - 1));
Irs_shaded = Isc / (exp(q*Voc/(Kb*A*Tk*Ns) - 1));
Io1 = exp(q*Ego*Tdif/(A*Kb*Tr*Tk)) * Irs * (Tk/Tr).^3;
Io1_shaded = exp(q*Ego*Tdif/(A*Kb*Tr*Tk)) * Irs_shaded * (Tk/Tr).^3;
%Io = Io1 * (exp(q*(Vpv+Ipv*Rs)/(Ns*A*Kb*Tk)) - 1);
%Ipv = Iph - Io - (Vpv+Ipv*Rs*Ns)/(Ns*Rp);

It = linspace(Iph_light, 0, 1000).';
Vt = zeros(size(It));
V1_last = 0; 
V2_last = 0;

for k = 1:numel(It)
    I = It(k);
    F1 = @(V) Iph_light - I - Io1 * (exp(q*(V + I*Rs)/(Ns*A*Kb*Tk)) - 1) - (V + I*Rs*Ns)/(Ns*Rp);
    V1 = fsolve(F1, V1_last);

    if I <= Iph_shaded
        F2 = @(V) Iph_shaded - I - Io1_shaded * (exp(q*(V + I*Rs)/(Ns*A*Kb*Tk)) - 1) - (V + I*Rs*Ns)/(Ns*Rp);
        V2 = fsolve(F2, V2_last);
    else
        V2 = -(Vf + I*Ron); 
    end
    Vt(k) = V1 + V2;  
    V1_last = V1;            
    V2_last = V2;
end
Pt = Vt .* It;
%% figure
figure
plot(Vt,It,'LineWidth',2);
xlim([0 85])
xlabel('Voltage (V)'); 
ylabel('Current (A)'); 
title('I-V curve of two PV cells in series for different irradiance');
grid on

figure
plot(Vt,Pt,'LineWidth',2);
xlim([0 85])
xlabel('Voltage (V)'); 
ylabel('Power (W)'); 
title('P–V curve of two PV cells in series for different irradiance');
grid on