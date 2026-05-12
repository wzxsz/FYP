%% Optimization Toolbox is required to run the emulator
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
G_light = 1000;   % Operating irradiance level for normal illumination PV cell(W/m2)
G_shaded = 500;   % Operating irradiance level for shaded PV cell(W/m2)
Tr = 298.15;% Reference cell temperature(K) of 25°C
Tk = 298.15;% Operating temperature(K) suppose as 25°C 
Tdif = Tk - Tr;
A = 1.2;   % Diode ideality factor
Ego = 1.15;% Band gap energy for the silicon semiconductor eV
Rs = 1;    % Series connected resistance
Rp = 1e3; % Parallel connected resistance 
Vf  = 0.6; % Forward voltage for diode
Ron = 1e-2;% Forward resistance for diode

%{
% Intermediate formulas
Iph_light = (Isc+Ki*Tdif) * G / Gr;
Iph_shaded = (Isc+Ki*Tdif) * G1 / Gr;
Irs = Isc / (exp(q*Voc/(Kb*A*Tk*Ns) - 1));
Irs_shaded = Isc / (exp(q*Voc/(Kb*A*Tk*Ns) - 1));
Io1 = exp(q*Ego*Tdif/(A*Kb*Tr*Tk)) * Irs * (Tk/Tr).^3;
Io1_shaded = exp(q*Ego*Tdif/(A*Kb*Tr*Tk)) * Irs_shaded * (Tk/Tr).^3;
Io = Io1 * (exp(q*(Vpv+Ipv*Rs)/(Ns*A*Kb*Tk)) - 1);
Ipv = Iph - Io - (Vpv+Ipv*Rs*Ns)/(Ns*Rp);
%}
Irs = Isc / (exp(q*Voc/(Kb*A*Tk*Ns) - 1));
Io1 = exp(q*Ego*Tdif/(A*Kb*Tr*Tk)) * Irs * (Tk/Tr).^3;
opts = optimoptions('fsolve','Display','off');

nRow = 4;% Number of rows in the PV array
nCol = 4;% Number of columns in the PV array
nModule = nRow * nCol;% Total number of PV modules

cloudH = 2;   % Height of moving cloud mask (rows)
cloudW = 2;   % Width of moving cloud mask (columns)
t_list = 0:6; % Time steps for cloud movement
nTime = numel(t_list);% Total number of time steps
% Cell arrays to store row-wise I-V and P-V curves at each time
V_all = cell(nTime,nRow);
I_all = cell(nTime,nRow);
P_all = cell(nTime,nRow);
% Matrices to store row-wise MPP values at each time
Vmp_all = zeros(nTime,nRow);
Imp_all = zeros(nTime,nRow);
Pmp_all = zeros(nTime,nRow);
% Vector to store total maximum power at each time step
Ptotal_all = zeros(nTime,1);
% Cell array to store irradiance maps at each time step
Gmap_all = cell(nTime,1);
%% Main simulation loop over time
for tt = 1:nTime
    t = t_list(tt);
    Gmap = G_light * ones(nRow, nCol);% Start with a fully illuminated 4x4 irradiance map
    % Define the top and bottom row position of the moving cloud
    topRow = 1 - cloudH + t;
    bottomRow = topRow + cloudH - 1;
    overlapRows = max(1, topRow):min(nRow, bottomRow);% Find the rows of the cloud that overlap with the PV array
    if ~isempty(overlapRows)% Apply shading to the left cloudW columns within overlapped rows
        Gmap(overlapRows, 1:cloudW) = G_shaded;
    end
    Gmap_all{tt} = Gmap;% Save the irradiance map of current time step
    
    % Solve each row independently, each row is treated as one series string
    for r = 1:nRow
        G_row = Gmap(r,:).';% Extract irradiance values for the current row (4 modules) 
        Iph_vec = (Isc + Ki*Tdif) .* G_row / Gr;% Compute photocurrent of each module in this row
        It = linspace(max(Iph_vec), 0, 1000).';% Sweep the row current from short-circuit region to zero
        Vt = zeros(size(It));% Preallocate total row voltage vector
        V_last = zeros(nCol,1);% Store previous voltage solutions as better initial guesses for fsolve

    % Loop through all current points
    for k = 1:numel(It)
        I = It(k);
        Vsum = 0;% Total row voltage at this current
        % Solve each module voltage in the row
        for m = 1:nCol
            if I <= Iph_vec(m)% If row current does not exceed module photocurrent, the module stays on its intrinsic I-V curve
                Fm = @(V) Iph_vec(m) - I - Io1 * (exp(q*(V + I*Rs)/(Ns*A*Kb*Tk)) - 1) - (V + I*Rs*Ns)/(Ns*Rp);% Formula used to solve
                Vm = fsolve(Fm, V_last(m), opts);% Solve implicit module voltage using previous solution as initial guess
            else% If row current exceeds photocurrent, bypass diode turns on
                Vm = -(Vf + I*Ron);% Formula of diode voltage drop
            end
            Vsum = Vsum + Vm;% Add module voltage to obtain total row voltage
            V_last(m) = Vm;% Update previous voltage guess for next current point
        end
        Vt(k) = Vsum;% Store total row voltage for this current point
    end

    Pt = Vt .* It;% Compute row power curve
    
    % Find maximum power point (MPP) of this row
    [Pmp_now, idx] = max(Pt);
    Vmp_now = Vt(idx);
    Imp_now = It(idx);
    
    % Store row-wise I-V and P-V curves
    V_all{tt, r} = Vt;
    I_all{tt, r} = It;
    P_all{tt, r} = Pt;

    % Store row-wise MPP values
    Vmp_all(tt, r) = Vmp_now;
    Imp_all(tt, r) = Imp_now;
    Pmp_all(tt, r) = Pmp_now;
    end
    % Total maximum power at current time:sum of the maximum power of all rows
    Ptotal_all(tt) = sum(Pmp_all(tt,:));

    % Print row-wise and total MPP results
    fprintf('t = %d s\n', t);
    for r = 1:nRow
        fprintf('  Row %d -> Vmp = %.3f V, Imp = %.3f A, Pmp = %.3f W\n', r, Vmp_all(tt,r), Imp_all(tt,r), Pmp_all(tt,r));
    end
    fprintf('  Total power (sum of row MPPs) = %.3f W\n\n', Ptotal_all(tt));
end
%%
% Plot row-wise I-V curves over time
for r = 1:nRow
    figure
    hold on
    % Plot I-V curves of the same row at different time steps
    for tt = 1:nTime
        plot(V_all{tt,r}, I_all{tt,r}, 'LineWidth', 2)
    end
    grid on
    xlabel('Voltage (V)')
    ylabel('Current (A)')
    title(['I-V curves of Row ' num2str(r) ' at different times'])
    legend('t=0 s','t=1 s','t=2 s','t=3 s','t=4 s','t=5 s','t=6 s','Location','best')
    xlim([0 180])
end

% Plot row-wise P-V curves and MPPs over time
for r = 1:nRow
    figure
    hold on
    % Plot P-V curves of the same row at different time steps
    for tt = 1:nTime
        plot(V_all{tt,r}, P_all{tt,r}, 'LineWidth', 2)
        plot(Vmp_all(tt,r), Pmp_all(tt,r), 'o', 'MarkerSize', 7, 'LineWidth', 1.5)
    end
    grid on
    xlabel('Voltage (V)')
    ylabel('Power (W)')
    title(['P-V curves of Row ' num2str(r) ' at different times'])
    legend('t=0','MPP t=0', 't=1','MPP t=1', 't=2','MPP t=2', 't=3','MPP t=3', 't=4','MPP t=4', 't=5','MPP t=5', 't=6','MPP t=6','Location','best')
    xlim([0 180])
end

% Plot total maximum power point value over time
figure
bar(t_list, Ptotal_all)
grid on
xlabel('Time (s)')
ylabel('Total Maximum Power (W)')
title('Total Maximum Power at Each Time')