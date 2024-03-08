% AER E 344 Spring 2024 Lab 06 Part 2 Analysis
% Section 3 Group 3

clear, clc, close all;

% Initialize arrays
pitotTube_inH2O = [];
hotWireV = [];

% Loop through files
for fileIndex = 20:4:60
    % Generate filename
    filename = sprintf('%d.txt', fileIndex);

    % Read the file
    data = readtable(filename, 'HeaderLines', 5, 'Delimiter', '\t', 'Format', '%s%f%f', 'ReadVariableNames', false);

    % Extract columns 3 and 4 from the data table
    column3 = data.Var2;
    column4 = data.Var3;

    % Calculate mean values and append to arrays
    meanColumn3 = mean(column3);
    meanColumn4 = mean(column4);
    pitotTube_inH2O = [pitotTube_inH2O, meanColumn4];
    hotWireV = [hotWireV, meanColumn3];


end

% Display pitotTube_inH2O and hotWireV arrays
disp('Mean Values for Pitot Tube Pressure [inH2O]:');
disp(pitotTube_inH2O);
disp('Mean Values for Hot Wire Voltage [V]:');
disp(hotWireV);

% Variables
rho_air = 1.195; % [kg / m^3]

%Estimating the pitot tube V when the motor is off
motorFreqPercentage = 0.2:0.04:0.6; % % of max motor frequency
P1 = polyfit(motorFreqPercentage, pitotTube_inH2O, 1);
pitotTube_inH2Oat0 = P1(2); % pitot tube V when the motor is turned off

%Calculate v_pt
q_pt = (pitotTube_inH2O - pitotTube_inH2Oat0) .* 248.84; % [Pa]
v_pt = sqrt(2 .* q_pt ./ rho_air); % [m/s]

%Display v_T array
disp('Flow velocity [m/s]:')
disp(v_pt)

%Calculate and display curve fitting
P3 = polyfit(hotWireV, v_pt, 4);
x = hotWireV(1):0.01:hotWireV(11);
curveFit = P3(1) * x.^4 + P3(2) * x.^3 + P3(3) * x.^2 + P3(4) * x + P3(5);
fprintf('4th order polynomial curve fitting formula: \n   %.4fx^4 + %.4fx^3 + %.4fx^2 + %.4fx + %.4f \n\n', P3(1), P3(2), P3(3), P3(4), P3(5))

%Plot Hot Wire Voltage vs. Air speed
plot(hotWireV, v_pt, 'o')
hold on
plot(x, curveFit)
xlabel('Hot Wire Voltage [V]')
ylabel('Air speed [m/s]')
legend('Experimental data', '4th order polynomial curve fitting', 'Location', 'SouthEast')
grid
