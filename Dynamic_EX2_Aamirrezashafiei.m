%.....Dynamic_EX2_Amirreza shafiei_402449105.....
clc;
clear;

%% Given parameters:
omega0 = 1.0;          % Initial rotor speed in p.u.
delta0 = deg2rad(30);  % Initial rotor angle in radians
Pm = 1.0;              % Mechanical power in p.u.
H = 5.0;               % Inertia constant in seconds
faultStartTime = 1.0;  % Time when fault occurs in seconds
faultClearTime = 1.1;  % Time when fault is cleared in seconds
postFaultPe = 1.2;     % Electrical power after fault is cleared in p.u.
faultDuration = 0.1;   % Duration of the fault in seconds
stepSize = 0.01;       % Time step for integration

%% Create time array:
time = 0:stepSize:1.3;

%% Initialize arrays to store results:
deltaValues = zeros(size(time));
omegaValues = zeros(size(time));
PeValues = zeros(size(time));
PaccValues = zeros(size(time));

%% Set initial values:
deltaValues(1:find(time == 1.0)) = delta0;  
omegaValues(1:find(time == 1.0)) = omega0;  
PeValues(1:find(time == 1.0)) = 1.0;        

%% Function to calculate Pe based on time:
calculatePe = @(t, faultStartTime, faultClearTime, postFaultPe) ...
    (t < faultStartTime) * 1.0 + ...
    (t >= faultStartTime && t < faultClearTime) * 0.4 + ...
    (t >= faultClearTime) * postFaultPe;

%% Numerical integration using the trapezoidal method:
for i = find(time == 1.0):numel(time)-1
    t = time(i);
    Pe = calculatePe(t, faultStartTime, faultClearTime, postFaultPe);
    deltaPrev = deltaValues(i);
    omegaPrev = omegaValues(i);
    Pacc = Pm - Pe;

    deltaNext = deltaPrev + (stepSize / 2) * (omegaPrev + omegaPrev);
    omegaNext = omegaPrev + (stepSize / (2 * H)) * Pacc;

    deltaValues(i+1) = deltaNext;
    omegaValues(i+1) = omegaNext;
    PeValues(i+1) = Pe;
    PaccValues(i+1) = Pacc;
end

%% Analyze results:
stable = true;       % Assume the system is initially stable
criticalDelta = -1;  % Initialize critical delta as -1

for i = 1:numel(time)
    if abs(deltaValues(i)) > pi/6  % Define a threshold for stability
        fprintf('Time: %.2f s - Unstable condition\n', time(i));
        fprintf('δ: %.2f degrees\n', rad2deg(deltaValues(i)));
        fprintf('ω: %.2f p.u.\n', omegaValues(i));
        fprintf('P_{acc}: %.2f p.u.\n\n');
        stable = false;                  
        criticalDelta = deltaValues(i);  
        break;                           
    end
end

if stable
    fprintf('The system is stable within the defined threshold.\n');
else
    fprintf('Critical Delta: %.2f degrees\n', rad2deg(criticalDelta));
end

% Create plots
figure;

%% Plot Pe in terms of time:
subplot(4, 1, 1);
plot(time, PeValues, 'b', 'LineWidth', 1.2);
title('Electrical Power (P_e) vs. Time');
xlabel('Time (s)');
ylabel('P_e (p.u.)');
grid on;
legend('P_e (p.u.)');

%% Plot omega in terms of time:
subplot(4, 1, 2);
plot(time, omegaValues, 'g', 'LineWidth', 1.2);
title('Rotor Speed (ω) vs. Time');
xlabel('Time (s)');
ylabel('ω (p.u.)');
grid on;
legend('ω (p.u.)');

%% Plot Pacc in terms of time:
subplot(4, 1, 3);
plot(time, PaccValues, 'r', 'LineWidth', 1.2);
hold on;
yline(0, '--', 'Mechanical Power Balance', 'Color', 'k', 'LineWidth', 1.2);
title('Mechanical Power Deficit (P_{acc}) vs. Time');
xlabel('Time (s)');
ylabel('P_{acc} (p.u.)');
grid on;
legend('P_{acc} (p.u.)', 'Location', 'SouthEast');

%% Plot delta in terms of time:
subplot(4, 1, 4);
plot(time, rad2deg(deltaValues), 'm', 'LineWidth', 1.2);
title('Rotor Angle (δ) vs. Time');
xlabel('Time (s)');
ylabel('δ (degrees)');
grid on;
legend('δ (degrees)');
