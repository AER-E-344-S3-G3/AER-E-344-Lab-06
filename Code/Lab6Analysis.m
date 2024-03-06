clear; clc; close all;

%% Knowns
figure_dir = "../Figures/";
data_dir = "./Data/";
% Viscosity calculated at 21.2°C
% https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601
% .html
mu_air = 18.18e-6; % [Pa*s]
rho_air = 1.225; % [kg/m^3]
V_inf = 19.4; % [m/s} 
K = 1.1;
chord = 0.101; % [m] 
d_taps = 0.002; % [m] distance between taps

%% Part 1 Airfoil Wake 
aoa = [0 4 6 8 10 12];
n_aoa = length(aoa);
wake_files = dir(data_dir + '*deg.csv');
range = cat(2,(2:17),(35: 50),(68: 81)); % Columns containing pressure values from raw data
n_taps = length(range); % number of taps
position = (1:n_taps); % position of taps

PA_idx = 82;
PE_idx = 83;
C_p = zeros(n_aoa, n_taps); % Initialize matrix for pressure data
U_p = zeros(n_aoa, n_taps); % Initialize matrix for velocity data
C_D = zeros(1,n_aoa);
U_inf = zeros(1,n_aoa);
for a = 1:n_aoa % Iterate through files
    temp = mean(readmatrix(data_dir+wake_files(a).name)); % Read file & find average of columns
    count = 1;
    PA = temp(PA_idx);
    PE = temp(PE_idx);
    q = K *(PA - PE);
    U_inf(a) = sqrt(2 * q / rho_air );
    for i = 1 : width(temp) % Iterate through data
        if ismember(i,range) % Take columns with pressure data
            C_p(a, count) = (temp(i)-PE) / q ; % Cp at each tap
            U_p(a,count) = sqrt((temp(i) - PE) * 2 / rho_air); % Velocity at each tap
            count = count + 1;

        end 
    end
    figure(a)
    plot(position, C_p(a,:))
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "C_p Distribution of Airfoil Wake at " + aoa(a) +"^{\circ} AOA";
    title(title_str);
    xlabel("Position");
    ylabel("C_p");
    grid on;
    % saveas(gcf, figure_dir + title_str + ".svg");

% midpoint rule for approx integral
    for i = 1:range-1
        f1 = U_p(a, i)/U_inf(a) * (1 - U_p(a, i)/U_inf(a));
        f2 = U_p(a, i+1)/U_inf(a) * (1 - U_p(a, i+1)/U_inf(a));
        C_D(a) = C_D(a) + d_taps*((f1 + f2)/2);
    end
    C_D(a) = C_D(a) * 2 / chord;
end

figure(a + 1)
plot(aoa, C_D)
    fontname("Times New Roman");
    fontsize(12, "points");
    title_str = "C_D Distribution of Airfoil vs. AOA";
    title(title_str);
    xlabel("Angle of Attack [^{\circ}]");
    ylabel("C_D");
    grid on;

%% Part 2 Anemometer Calibration
%TODO