
%% Amirreza Shafiei/ ID: 402449105/ Date: 4030909

clc;
clear;
close all;

%% ---------------------Input Variables-------------------
Xd  = 1;
Xq  = 0.8;
Xpd = 0.1;
H   = 5;
Tdo = 8;
vt  = 1;
P0  = 0.8;
Q0  = 0.6;
F   = 50;
Ka  = 50;   % Exciter gain
Ta  = 0.4;  % Exciter time


%% ---------------------Time and Value Control-------------------
t_fault  = 2;         % Fault start time
t_clear  = 2.1;       % Fault clearance time
Sim_time = 50;        % Total simulation time
t_max    = Sim_time;  % Maximum simulation time
dt       = 0.02;

ICCT = 1;  % Flag to activate (1) or deactivate (0) CCT calculation

t_ref    = 200;       
Vref2    = 1.3;

t_m      = 300;      
Pm2      = 0.7;

%% --------------------- Initialization --------------------
w0        = 1;
S         = P0 + 1j * Q0;
I0        = (conj(S)) / (conj(vt));
oa        = vt + 1j * Xq * I0;
phi       = angle(oa);
vq0       = abs(vt) * cos(phi);
Id_0      = abs(I0) * sin(-angle(I0) + phi);
Eq_0      = vq0 + Xpd * Id_0;
Xe_PreF   = 0.4; 
Xe        = 0.4;
Xlines    = Xe / 2;
Xt        = 0.2;
XT        = Xt + Xlines;
V_infinit = vt - 1j * Xe_PreF * I0;
Vinf      = abs(V_infinit);
Delta_0   = (phi + abs(angle(V_infinit))) * (180) / pi;
Pm        = P0;
Efd_0     = (Xd - Xpd) * Id_0 + Eq_0;
Iq_0      = abs(I0) * cos(abs(angle(I0)) + phi);
Vd_0      = abs(vt) * cos((pi / 2) - phi);
P0        = Vd_0 * Id_0 + vq0 * Iq_0;
It_0      = abs(I0);
Vref      = abs(vt) + Efd_0 / Ka;
Vq_0      = abs(vt) * cos(phi);

%% --------------------- Simulation Phases ------------------------

%============================ Pre-Fault ==========================
t_PreF = 0:dt:t_fault-dt;
w(1:length(t_PreF)) = w0;
Delta(1:length(t_PreF)) = Delta_0;
Eq(1:length(t_PreF)) = Eq_0;
Efd(1:length(t_PreF)) = Efd_0;
Pe(1:length(t_PreF)) = P0;
Vq(1:length(t_PreF)) = Vq_0;
Vd(1:length(t_PreF)) = Vd_0;
Vt(1:length(t_PreF)) = vt;
Id(1:length(t_PreF)) = Id_0;
Iq(1:length(t_PreF)) = Iq_0;
I(1:length(t_PreF)) = It_0;

%============================ During Fault ==========================
t_DF = t_fault:dt:t_clear-dt;
L = length(t_PreF);
Dt = dt / 2;
S1 = dt / (2 * Ta);
S2 = dt / (2 * Tdo);
S3 = 1 + S1;
S4 = 1 - S1;
S5 = 1 + S2;
S6 = 1 - S2;
S7 = pi / 180;

for i = 1:length(t_DF)
    if t_DF(i) >= t_ref && t_DF(i) <= t_max
        Vref = Vref2;
    end
    if t_DF(i) >= t_m && t_DF(i) <= t_max
        Pm = Pm2;
    end

    Pe(L + i) = 0;
    w(L + i) = w(L + i - 1) + Dt * (w0 / (2 * H)) * (2 * Pm - (2 * Pm - Pe(L + i - 1) - Pe(L + i - 1)));
    Delta(L + i) = Delta(L + i - 1) + Dt * (w(L + i) + w(L + i - 1) - 2 * w0) * 18000;
    Efd(L + i) = (Efd(L + i - 1) * S4 + S1 * (2 * Ka * Vref)) / S3;
    Eq(L + i) = (Eq(L + i - 1) * S6 + S2 * (Efd(L + i - 1) + Efd(L + i))) / S5;
    Vq(L + i) = 0;
    Vd(L + i) = 0;
    Vt(L + i) = 0;
    Id(L + i) = 0;
    Iq(L + i) = 0;
    I(L + i) = 0;
end

%============================ Post Fault ==========================
t_PF = t_clear:dt:Sim_time;
j = length(t_DF);

for i = 1:length(t_PF)
    if t_PF(i) >= t_ref && t_PF(i) <= t_max
        Vref = Vref2;
    end
    if t_PF(i) >= t_m && t_PF(i) <= t_max
        Pm = Pm2;
    end
%------------------------- First Iteration ---------------------------
    w1 = w(L + j + i - 1) + Dt * (w0 / (2 * H)) * (2 * Pm - Pe(L + j + i - 1) - Pe(L + j + i - 1));
    Delta1 = Delta(L + j + i - 1) + Dt * (w(L + j + i - 1) + w(L + j + i - 1) - 2 * w0) * 18000;
    Efd1 = (Efd(L + j + i - 1) * S4 + S1 * (2 * Ka * Vref - Ka * (Vt(i) + Vt(i)))) / S3;
    Eq1 = (Eq(L + j + i - 1) * S6 + S2 * (Efd(L + j + i - 1) + Efd(L + j + i - 1) - (Xd - Xpd) * (Id(i) + Id(i)))) / S5;
    Vq1 = Vinf * cos(Delta1 * S7) + (XT / (Xpd + XT)) * (Eq1 - Vinf * cos(Delta1 * S7));
    Vd1 = (Xq / (Xq + XT)) * Vinf * sin(Delta1 * S7);
    Vt1 = sqrt(Vq1^2 + Vd1^2);
    Id1 = (Eq1 - Vinf * cos(Delta1 * S7)) / (Xpd + XT);
    Iq1 = Vinf * sin(Delta1 * S7) / (Xq + XT);
    Pe1 = Vd1 * Id1 + Vq1 * Iq1;

%------------------------- Second Iteration ---------------------------

    w(L + j + i) = w(L + j + i - 1) + Dt * (w0 / (2 * H)) * (2 * Pm - Pe1 - Pe(L + j + i - 1));
    Delta(L + j + i) = Delta(L + j + i - 1) + Dt * (w1 + w(L + j + i - 1) - 2 * w0) * 18000;
    Efd(L + j + i) = (Efd(L + j + i - 1) * S4 + S1 * (2 * Ka * Vref - Ka * (Vt1 + Vt(i)))) / S3;
    Eq(L + j + i) = (Eq(L + j + i - 1) * S6 + S2 * (Efd1 + Efd(L + j + i - 1) - (Xd - Xpd) * (Id1 + Id(i)))) / S5;

    Vq(i + 1) = Vinf * cos(Delta(L + j + i) * S7) + (XT / (Xpd + XT)) * (Eq(L + j + i) - Vinf * cos(Delta(L + j + i) * S7));
    Vd(i + 1) = (Xq / (Xq + XT)) * Vinf * sin(Delta(L + j + i) * S7);
    Vt(i + 1) = sqrt(Vq(i + 1)^2 + Vd(i + 1)^2);
    Id(i + 1) = (Eq(L + j + i) - Vinf * cos(Delta(L + j + i) * S7)) / (Xpd + XT);
    Iq(i + 1) = Vinf * sin(Delta(L + j + i) * S7) / (Xq + XT);
    Pe(L + j + i) = Vd(i + 1) * Id(i + 1) + Vq(i + 1) * Iq(i + 1);
    
end

  %% --------------------- Critical Clearing Time (CCT) Calculation ------------------------
if ICCT == 1
% Initialize CCT Calculation

t_clear = t_fault + 0.134; % Initial guess
deltat_DF = 0.134; % Fault clearing time increment
max_iterations = 100;
iteration = 0;
if t_PF(i) >= t_ref
   t_clear = t_fault + 0.113;
end
    if t_PF(i) >= t_m
       t_clear = t_fault + 0.099;
    end
while iteration < max_iterations

    iteration = iteration + 1;

%============================ Pre-Fault ==========================
    t_PreF = 0:dt:t_fault-dt;
    w(1:length(t_PreF)) = w0;
    Delta(1:length(t_PreF)) = Delta_0;
    Eq(1:length(t_PreF)) = Eq_0;
    Efd(1:length(t_PreF)) = Efd_0;
    Pe(1:length(t_PreF)) = P0;
    Vq(1:length(t_PreF)) = Vq_0;
    Vd(1:length(t_PreF)) = Vd_0;
    Id(1:length(t_PreF)) = Id_0;
    Iq(1:length(t_PreF)) = Iq_0;
    I(1:length(t_PreF)) = It_0;

%============================ During Fault ==========================
    for i = 1:length(t_DF)
        if t_PF(i) >= t_ref
        Vref = 1.3;
        end
        if t_PF(i) >= t_m
        Pm = 0.7;
        end
        Pe(L + i) = 0;
        w(L + i) = w(L + i - 1) + Dt * (w0 / (2 * H)) * (2 * Pm - (2 * Pm - Pe(L + i - 1) - Pe(L + i - 1)));
        Delta(L + i) = Delta(L + i - 1) + Dt * (w(L + i) + w(L + i - 1) - 2 * w0) * 18000;
        Efd(L + i) = (Efd(L + i - 1) * S4 + S1 * (2 * Ka * Vref)) / S3;
        Eq(L + i) = (Eq(L + i - 1) * S6 + S2 * (Efd(L + i - 1) + Efd(L + i))) / S5;
        Vq(L + i) = 0;
        Vd(L + i) = 0;
        Id(L + i) = 0;
        Iq(L + i) = 0;
        I(L + i) = 0;
    end

%============================ Post Fault ==========================
    stable = true;
    for i = 1:length(t_PF)
        if t_PF(i) >= t_ref
        Vref = 1.3;
        end
        if t_PF(i) >= t_m
        Pm = 0.7;
        end
        %------------------------- First Iteration ---------------------------

        w1 = w(L + j + i - 1) + Dt * (w0 / (2 * H)) * (2 * Pm - Pe(L + j + i - 1) - Pe(L + j + i - 1));
        Delta1 = Delta(L + j + i - 1) + Dt * (w(L + j + i - 1) + w(L + j + i - 1) - 2 * w0) * 18000;
        Efd1 = (Efd(L + j + i - 1) * S4 + S1 * (2 * Ka * Vref - Ka * (Vt(i) + Vt(i)))) / S3;
        Eq1 = (Eq(L + j + i - 1) * S6 + S2 * (Efd(L + j + i - 1) + Efd(L + j + i - 1) - (Xd - Xpd) * (Id(i) + Id(i)))) / S5;
        Vq1 = Vinf * cos(Delta1 * S7) + (XT / (Xpd + XT)) * (Eq1 - Vinf * cos(Delta1 * S7));
        Vd1 = (Xq / (Xq + XT)) * Vinf * sin(Delta1 * S7);
        Id1 = (Eq1 - Vinf * cos(Delta1 * S7)) / (Xpd + XT);
        Iq1 = Vinf * sin(Delta1 * S7) / (Xq + XT);
        Pe1 = Vd1 * Id1 + Vq1 * Iq1;

       %------------------------- Second Iteration ---------------------------

        w(L + j + i) = w(L + j + i - 1) + Dt * (w0 / (2 * H)) * (2 * Pm - Pe1 - Pe(L + j + i - 1));
        Delta(L + j + i) = Delta(L + j + i - 1) + Dt * (w1 + w(L + j + i - 1) - 2 * w0) * 18000;
        Efd(L + j + i) = (Efd(L + j + i - 1) * S4 + S1 * (2 * Ka * Vref - Ka * (Vt1 + Vt(i)))) / S3;
        Eq(L + j + i) = (Eq(L + j + i - 1) * S6 + S2 * (Efd1 + Efd(L + j + i - 1) - (Xd - Xpd) * (Id1 + Id(i)))) / S5;
        Vq(i + 1) = Vinf * cos(Delta(L + j + i) * S7) + (XT / (Xpd + XT)) * (Eq(L + j + i) - Vinf * cos(Delta(L + j + i) * S7));
        Vd(i + 1) = (Xq / (Xq + XT)) * Vinf * sin(Delta(L + j + i) * S7);
        Id(i + 1) = (Eq(L + j + i) - Vinf * cos(Delta(L + j + i) * S7)) / (Xpd + XT);
        Iq(i + 1) = Vinf * sin(Delta(L + j + i) * S7) / (Xq + XT);
        Pe(L + j + i) = Vd(i + 1) * Id(i + 1) + Vq(i + 1) * Iq(i + 1);
        if max(Delta(L+1:end)) < 180
           t_clearCCT = t_clear + deltat_DF;
            stable = false;
            break;
        end
    end

    if stable
       t_clearCCT = t_clear + deltat_DF; % Increase fault clearing time
    else
        t_clearCCT = t_clear - deltat_DF; % Step back to stable time
        CCT = (t_clear - t_fault) * 1000; % Convert to milliseconds
        fprintf('Critical Clearing Time (CCT): %.2f ms\n', CCT);
        break;
    end
end

if iteration == max_iterations
    disp('CCT calculation did not converge.');
end

else
    fprintf('CCT calculation is inactive.\n');
end


%% --------------------- Save Results to Excel ------------------------
time = 0:dt:Sim_time;
min_length = min([length(time), length(Delta), length(w), length(Pe), length(Eq), length(Efd), length(Vt)]);
output_filename = 'Simulation_Results.xlsx';

time_trimmed = time(1:min_length);
Delta_trimmed = Delta(1:min_length);
w_trimmed = w(1:min_length);
Pe_trimmed = Pe(1:min_length);
Eq_trimmed = Eq(1:min_length);
Efd_trimmed = Efd(1:min_length);
Vt_trimmed = Vt(1:min_length);

input_table = table({'Xd'; 'Xq'; 'Xpd'; 'H'; 'Tdo'; 'vt'; 'P0'; 'Q0'; 'Ka'; 'Ta'; 'dt'; ...
                     't_fault'; 't_clear'; 'Sim_time'; 't_ref'; 't_m'}, ...
                    [Xd; Xq; Xpd; H; Tdo; vt; P0; Q0; Ka; Ta; dt; ...
                     t_fault; t_clear; Sim_time; t_ref; t_m], ...
                    'VariableNames', {'Parameter', 'Value'});

results_table = table(time_trimmed', Delta_trimmed', w_trimmed', Pe_trimmed', Eq_trimmed', Efd_trimmed', Vt_trimmed', ...
    'VariableNames', {'Time', 'Delta', 'RotorSpeed', 'ElectricalPower', 'Eq', 'Efd', 'Vt'});

if ICCT == 1
    results_table.CCT = nan(height(results_table), 1); % Add a CCT column filled with NaN
    results_table.CCT(1) = CCT; % Set the first row's CCT value
end

writetable(input_table, output_filename, 'Sheet', 'Input Parameters');

writetable(results_table, output_filename, 'Sheet', 'Output Results');





%% --------------------- Display Results ------------------------
figure;
sgtitle(sprintf('System Parameters:\nXd = %.2f, Xq = %.2f, Xpd = %.2f, H = %.2f, Tdo = %.2f, vt = %.2f, P0 = %.2f, Q0 = %.2f, Ka = %.2f, Ta = %.2f', ...
                Xd, Xq, Xpd, H, Tdo, vt, P0, Q0, Ka, Ta));


subplot(3, 2, 1);
plot(time, Delta, 'LineWidth', 1.5);
xlabel('\fontname{times new roman}Time (s)','Fontsize',15);
ylabel('\fontname{times new roman}\delta','Fontsize',15);
title('Rotor Angle vs. Time');

subplot(3, 2, 2);
plot(time, w, 'LineWidth', 1.5);
xlabel('\fontname{times new roman}Time (s)','Fontsize',15);
ylabel('\fontname{times new roman}W','Fontsize',15);
title('Rotor Speed (\omega) vs. Time');


subplot(3, 2, 3);
plot(time, Pe, 'LineWidth', 1.5);
xlabel('\fontname{times new roman}Time (s)','Fontsize',15);
ylabel('\fontname{times new roman}Pe','Fontsize',15);
title('Electrical Power vs. Time');

subplot(3, 2, 4);
plot(time, Eq, 'LineWidth', 1.5);
xlabel('\fontname{times new roman}Time (s)','Fontsize',15);
ylabel('\fontname{times new roman}E^\prime_q','Fontsize',15);
title('Eq vs. Time');

subplot(3, 2, 5);
plot(time, Efd, 'LineWidth', 1.5);
xlabel('\fontname{times new roman}Time (s)','Fontsize',15);
ylabel('\fontname{times new roman}E_{fd}','Fontsize',15);
title('Field Voltage vs. Time');

subplot(3, 2, 6);
plot(time_trimmed, Vt_trimmed, 'LineWidth', 1.5);
xlabel('\fontname{times new roman}Time (s)','Fontsize',15);
ylabel('\fontname{times new roman}V_t');
title('Terminal Voltage (V_t) vs. Time');



