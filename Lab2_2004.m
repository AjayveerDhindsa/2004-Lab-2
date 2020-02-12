%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Convert 2D to 3D data for coefficient of lift and drag.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% House Keeping
clearvars
close all
clc

%Defining Variables
a_Attack = -5:12; %Entire Plane
C_L = [-0.32438 -0.21503 -0.10081 0.010503 0.12155 0.24163 0.34336 0.45256 0.56037 0.66625 0.76942 0.86923 0.96386 1.0441 1.0743 1.0807 1.0379 1.034]; %Entire Plane
C_D = [0.044251 0.033783 0.028627 0.025864 0.024643 0.025099 0.025635 0.02766 0.030677 0.034855 0.040403 0.04759 0.057108 0.070132 0.090921 0.11193 0.13254 0.15645]; %Entire Plane

Re30data  = readcell('AG18_T1_Re0.030_M0.00_N9.0.txt');
Re40data  = readcell('AG18_T1_Re0.040_M0.00_N9.0.txt');
Re60data  = readcell('AG18_T1_Re0.060_M0.00_N9.0.txt');
Re80data  = readcell('AG18_T1_Re0.080_M0.00_N9.0.txt');
Re100data = readcell('AG18_T1_Re0.100_M0.00_N9.0.txt');

Data = {Re30data, Re40data, Re60data, Re80data, Re100data};

Data_Alpha = {Data{1}(8:end,1), Data{2}(8:end,1), Data{3}(8:end,1), Data{4}(8:end,1), Data{5}(8:end,1)};

Data_Cl = {Data{1}(8:end,2), Data{2}(8:end,2), Data{3}(8:end,2), Data{4}(8:end,2), Data{5}(8:end,2)};

Data_Cd = {Data{1}(8:end,3), Data{2}(8:end,3), Data{3}(8:end,3), Data{4}(8:end,3), Data{5}(8:end,3)};

%Plotting angle of attack and coefficients for lift and drag

%2D Plots
figure(1);
hold on;
plot(cell2mat(Data_Alpha{1}(:)), cell2mat(Data_Cl{1}(:)))
plot(cell2mat(Data_Alpha{1}(:)), cell2mat(Data_Cd{1}(:)))
title('Angle of Attack vs. Lift and Drag Coefficients of AG18 Airfoil')

% Determining the slope of each curve
% Using y=mx+b to find a_0 on 2D wing plot
a_0 = polyfit(cell2mat(Data_Alpha{1}(:)), cell2mat(Data_Cl{1}(:)), 1);
a_02D = a_0(1);

%Using Equation 3.1 to solve for a for the wing
e = 0.54; % Calculated Value
AR = 14.4;
a_wing3D = a_02D/(1+((57.3*a_02D)/(pi*e*AR)));

%Using y=mx+b to find a_0 on 2D airfoil plot
A = polyfit(a_Attack(1:12),C_L(1:12),1);
a_0airfoil = A(1);

%Using Equation 3.1 to solve for a_0 on 2D airfoil plot
a_airfoil3D = a_0airfoil/(1+(57.3*a_0airfoil)/(pi*e*AR)); %3D slope

%Calculating the angle of attack when lift = 0
%Using y=mx+b to solve for a_attack3D for airfoil
%Lift = 1/2 *C_l*rho*V^2*S
a_attack3Dwing = -cell2mat(Data_Cl{1}(67))/a_wing3D;

%Coefficient of Lift for 3D finite wing
CL_wing3D = a_wing3D.*(cell2mat(Data_Alpha{1}(:))-a_attack3Dwing);

%Coefficient of Drag for 3D finite wing
CD_wing3D = cell2mat(Data_Cd{1}(:))+(CL_wing3D.^2/(pi*e*AR));

%% Question 1

figure(2);
hold on;
plot(cell2mat(Data_Alpha{1}(:)),CL_wing3D)
% Overlaying the plots
plot(cell2mat(Data_Alpha{1}(:)),cell2mat(Data_Cl{1}(:)))

legend('3D Finite Wing','AG18 Airfoil','Location','SouthEast')
title('Angle of Attack vs. Lift of Finite 3D Wing, AG18 Airfoil')

%% Question 3

L_D = CL_wing3D./CD_wing3D;
CL_CD = C_L./C_D;
% Plotting L/D vs aoa of wing and tempest data
figure(3);
plot(cell2mat(Data_Alpha{1}(:)), L_D)
title('Angle of Attack vs. Lift Over Drag of 3D Finite Wing')
xlabel('Angle of Attack [degrees]')
ylabel('L/D')
% Overlaying Tempest UAS Data for Question 3

% Estimated Max L/D and V/AoA at Max L/D
for i = 1:length(L_D)
    if L_D(i) == max(L_D)
        Max_LD = max(L_D);
        Max_AoA = cell2mat(Data_Alpha{1}(i));
        break
        % Max velocity found at end with assumption of flying at 1.8km height 
    end
end
h = 7; % Height of Garage = 7m
Max_Range = Max_LD*h;
%%
% Using equation 3.6 to solve for CD_min to eventually solve 3.4a/b for CD
C_fe = 0.006; % 0.027 / Re^(1/7)
S_wet = (pi*0.015*0.2) + (2*0.06*0.4) + (0.05*0.03);
S_ref = 0.024;
CD_min = (C_fe * S_wet) / (S_ref);

% Solving 3.3b for k1 to solve 3.4a/b
k1 = 1 / (pi*e*AR);
% Initial Estimation of Oswalds efficiency to solve for 3.4a/b
e_0 = 1 / (1.05 + 0.007*pi*S_ref);

% Solving 3.3b for k1 to solve 3.4a/b
k1 = 1 / (pi*e_0*AR);

% Solving 3.5 for CL_minD for 3.4a/b eventually
CL_minD = a_wing3D*(CD_wing3D(5) - a_attack3Dwing); 

% Solving 3.4a for CD of entire aircraft
CD_Aircraft = CD_min + k1*(CL_wing3D - CL_minD).^2;
CD_0 = CD_min + k1*CL_minD.^2;

%Calculating The coefficients of Drag
k2 = -2*k1*CL_minD;
CD = CD_0 + k1*CL_wing3D.^2 + k2*CL_wing3D;

%Plotting new drag polar
figure(4)
hold on
plot(cell2mat(Data_Alpha{1}(:)),CD)
title('Drag Coefficient for the Entire Aircraft')
xlabel('Angle of Attack [degrees]')
ylabel('Drag Coefficient')
hold off

%% Weight Calculations
rho_wood_1 = 600; % [kg/m^3]
rho_wood_2 = 160; % [kg/m^3]
rho_wing = 50; % [kg/m^3]
W_payload = 0.160; % [kg]

W_fuselage = rho_wood_2 * (pi*(0.03^2)*0.35);
W_horztail = rho_wood_2 * (0.05*0.1*0.005);
W_verttail = rho_wood_2 * (0.05*0.05*0.005);
W_wing = rho_wing * (2*0.06*0.4*0.005);
W_wingrod_horz = rho_wood_2 * (pi*(0.01^2)*0.4*2*2);
W_wingrod_vert = rho_wood_2 * (pi*(0.01^2)*0.06*5*2);
W_ballast = ((W_fuselage*0.005) + ((W_verttail+W_horztail)*0.175))/(0.17);

W_total = W_fuselage + W_horztail + W_verttail + W_wing + W_wingrod_horz + W_wingrod_vert + W_payload + W_ballast; % [kg]

%%
% Now that we have CD_0 we can do max range and max endurance calculations
rho_air = 1.058; %[kg/m^3]

%determing angle of attack for zero lift for 3D airfoil
eq = polyfit(cell2mat(Data_Alpha{1}(:)), CL_wing3D,1);
a_L0 = eq(1);
b = eq(2);

% Max Range Finite Wing
CL_max = sqrt(CD_0/k1);
V_max = sqrt(2*W_total./(rho_air.*S_ref.*CL_max)); %[m/s]
a_attackL0 = (CL_max - b)/a_L0;
a_attackVmax = (CL_max/a_airfoil3D) - a_attackL0;

% Max Endurance Finite Wing
CL_max_E = sqrt(3*CD_0/k1);
V_max_E = sqrt(2*W_total./(rho_air.*S_ref.*CL_max_E)); %[m/s]
a_attackL0_E = (CL_max_E - b)/a_L0;
a_attackVmax_E = (CL_max_E/a_airfoil3D) - a_attackL0_E;

% V at Max L/D
V_LD = sqrt(2*W_total./(rho_air.*S_ref.*CL_wing3D(i)));

%% Stability
% Horizontal Tail Volume Coefficient (4.3.1)
V_H = ((0.04*0.1)*(0.18))/(S_ref*0.06);
% ~0.4667, Within 0.3 - 0.6 range

% Lateral Directional Stability
V_V = ((0.06*0.05)*(0.175))/(S_ref*0.83);
% ~0.031, Within 0.02 - 0.05 range

% Spiral Parameter
B = (0.18*6.5)/(0.83*CL_wing3D(i));
% 5.2, Greater than 5 so stable spiral