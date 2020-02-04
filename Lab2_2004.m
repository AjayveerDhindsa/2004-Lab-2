%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Purpose: Convert 2D to 3D data for coefficient of lift and drag.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% House Keeping
clear all
close all
clc

%Defining Variables
a_Attack = -5:12; %Entire Plane
C_L = [-0.32438 -0.21503 -0.10081 0.010503 0.12155 0.24163 0.34336 0.45256 0.56037 0.66625 0.76942 0.86923 0.96386 1.0441 1.0743 1.0807 1.0379 1.034]; %Entire Plane
C_D = [0.044251 0.033783 0.028627 0.025864 0.024643 0.025099 0.025635 0.02766 0.030677 0.034855 0.040403 0.04759 0.057108 0.070132 0.090921 0.11193 0.13254 0.15645]; %Entire Plane
a_attack = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15]; %Just Wings

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
%3D Plots
figure(1);
hold on;
plot(a_Attack,C_L)
plot(a_Attack,C_D)
title('Angle of Attack vs. Lift and Drag Coefficients of Tempest UAS')
xlabel('Angle of Attack [degrees]')
ylabel('Coefficents for Lift and Drag')

%2D Plots
figure(2);
hold on;
plot(cell2mat(Data_Alpha{3}(:)), cell2mat(Data_Cl{3}(:)))
plot(cell2mat(Data_Alpha{3}(:)), cell2mat(Data_Cd{3}(:)))
title('Angle of Attack vs. Lift and Drag Coefficients of Airfoil')
xlabel('Angle of Attack [degrees]')
ylabel('Coefficents for Lift and Drag')

% Determining the slope of each curve
% Using y=mx+b to find a_0 on 2D wing plot
a_0 = polyfit(cell2mat(Data_Alpha{3}(:)), cell2mat(Data_Cl{3}(:)), 1);
a_02D = a_0(1);

%Using Equation 3.1 to solve for a for the wing
e = 0.9387; %Given Value
AR= 16.5; %
a_wing3D = a_02D/(1+((57.3*a_02D)/(pi*e*AR)));

%Sparcing data for 2D wing plot
A_data = a_attack(1:12);
C_Ldata = C_L(1:12);
C_Ddata = C_D(1:12);

%Using y=mx+b to find a_0 on 2D airfoil plot
A = polyfit(A_data,C_Ldata,1);
a_0airfoil = A(1);

%Using Equation 3.1 to solve for a_0 on 2D airfoil plot
a_airfoil3D = a_0airfoil/(1+(57.3*a_0airfoil)/(pi*e*AR)); %3D slope

%Calculating the angle of attack when lift = 0
%Using y=mx+b to solve for a_attack3D for airfoil
%Lift = 1/2 *C_l*rho*V^2*S
a_attack3Dwing = -C_l(6)/a_wing3D;

%Coefficient of Lift for 3D finite wing
CL_wing3D = a_wing3D.*(a_attack-a_attack3Dwing);

%Coefficient of Drag for 3D finite wing
CD_wing3D = C_d+(CL_wing3D.^2/(pi*e*AR));

%% Question 1

figure(3);
hold on;
plot(a_attack,CL_wing3D)
% plot(a_attack,CD_wing3D)

% Overlaying the 
plot(a_attack,C_l)
plot(a_Attack,C_L)

legend('3D Finite Wing', 'MH 32 Airfoil', 'Tempest UAS')
%title('Angle of Attack vs. Lift and Drag Coefficients of 3D Finite Wing')
title('Angle of Attack vs. Lift of Finite 3D Wing, Mh 32 Airfoil, and Tempest UAS')
xlabel('Angle of Attack [degrees]')
ylabel('Coefficent of Lift')

%% Question 3

L_D = CL_wing3D./CD_wing3D;
CL_CD = C_L./C_D;
% Plotting L/D vs aoa of wing and tempest data
figure(5);
hold on;
plot(a_attack, L_D)
title('Angle of Attack vs. Lift Over Drag of 3D Finite Wing and Tempest UAS')
xlabel('Angle of Attack [degrees]')
ylabel('L/D')
% Overlaying Tempest UAS Data for Question 3
plot(a_Attack, CL_CD)
legend('3D Finite Wing', 'Tempest UAS')
hold off

% Calculating percent difference between two curves
CL_Difference = (((CL_wing3D(1:18) - C_L) ./ C_L) * 100);
CL_Difference_Avg = mean(CL_Difference);
% Max percent difference is at 4 degrees which is 725% difference, avg of
% 37%

% Estimated Max L/D and V/AoA at Max L/D
for i = 1:length(L_D)
    if L_D(i) == max(L_D)
        Max_LD = max(L_D);
        Max_AoA = i - 6;
        % Max velocity found at end with assumption of flying at 1.8km height 
    end
end

for i = 1:length(CL_CD)
    if CL_CD(i) == max(CL_CD)
        Max_LD_UAS = max(CL_CD);
        Max_AoA_UAS = i - 6;
        % Max velocity found at end with assumption of flying at 1.8km height 
    end
end
%%
% Using equation 3.6 to solve for CD_min to eventually solve 3.4a/b for CD
C_fe = 0.004;
S_wet = 2.327;
S_ref = 0.698;
CD_min = (C_fe * S_wet) / (S_ref);

% Initial Estimation of Oswalds efficiency to solve for 3.4a/b
s = 0.995;
S = 0.698;
e_0 = 1 / (1.05 + 0.007*pi*S);

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
figure(6)
hold on
plot(a_attack,CD)
title('Drag Coefiificent for the Entire Aircraft')
xlabel('Angle of Attack [degrees]')
ylabel('Draf Coefficient')
hold off

% Now that we have CD_0 we can do max range and max endurance calculations
height = 1800; %[m]
GTOW = 6.4*9.81; %[kg]
rho = 1.026937; %[kg/m^3]
S = 0.698; %[m^2] @ 1.8 km

%determing angle of attack for zero lift for 3D airfoil
eq = polyfit(a_attack, CL_wing3D,1);
a_L0 = eq(1);
b = eq(2);

% Max Range Finite Wing
CL_max = sqrt(CD_0/k1);
V_max = sqrt(2*GTOW./(rho.*S.*CL_max)); %[m/s]
a_attackL0 = (CL_max - b)/a_L0;
a_attackVmax = (CL_max/a_airfoil3D) - a_attackL0;

% Max Endurance Finite Wing
CL_max_E = sqrt(3*CD_0/k1);
V_max_E = sqrt(2*GTOW./(rho.*S.*CL_max_E)); %[m/s]
a_attackL0_E = (CL_max_E - b)/a_L0;
a_attackVmax_E = (CL_max_E/a_airfoil3D) - a_attackL0_E;

%determing angle of attack for zero lift for 3D airfoil
eq = polyfit(a_Attack, C_L,1);
a_L0_UAS = eq(1);
b_UAS = eq(2);

% Max Range of Tempest
V_max_UAS = sqrt(2*GTOW./(rho.*S.*C_L(10))); %[m/s]
a_attackL0_UAS = (C_L(10) - b_UAS)/a_L0_UAS;
a_attackVmax_UAS = (C_L(10)/a_L0_UAS) - a_attackL0_UAS;

% Max Endurance of Tempest
C_D_E = k1*C_L(10);
C_L_E = sqrt((3*C_D_E)/k1);
V_max_UAS_E = sqrt(2*GTOW./(rho.*S.*C_L_E)); %[m/s]
a_attackL0_UAS_E = (C_L_E - b_UAS)/a_L0_UAS;
a_attackVmax_UAS_E = (C_L_E/a_L0_UAS) - a_attackL0_UAS;

% Percent Difference between 3D Finite Wing and Tempest
Finite_Tempest_Difference = 100*abs((V_max - V_max_UAS)/V_max_UAS);
Finite_Tempest_Difference_AoA = 100*abs((a_attackVmax - a_attackVmax_UAS)/a_attackVmax_UAS);
Finite_Tempest_Difference_E = 100*abs((V_max_E - V_max_UAS_E)/V_max_UAS_E);
Finite_Tempest_Difference_AoA_E = 100*abs((a_attackVmax_E - a_attackVmax_UAS_E)/a_attackVmax_UAS_E);

% V at Max L/D
V_LD = sqrt(2*GTOW./(rho.*S.*CL_wing3D(7)));
V_LD_UAS = sqrt(2*GTOW./(rho.*S.*C_L(10)));
%% Question 2
figure(7);
hold on;
%plot(a_attack,CD_wing3D)

% Overlaying the 
plot(C_l,C_d)
plot(C_L,C_D)

legend('MH 32 Airfoil', 'Tempest UAS')
title('Drag Polars of Mh 32 Airfoil and Tempest UAS')
xlabel('Coefficent of Lift')
ylabel('Coefficent of Drag')

figure(8);
hold on;
plot(CL_wing3D,CD_wing3D)

% Overlaying the plots

plot(C_L,C_D)
plot(CL_wing3D,CD)


title('Comparasons of Drag Polar')
xlabel('Coefficent of Lift')
ylabel('Coefficent of Drag')
legend('Finite Wing Drag Polar', 'Tempest UAS Data','Whole Aircraft Drag Poalr')
hold off

figure(9);
hold on;

comp = CD(1:18);
differ = C_D - comp;
differ = (differ./C_D).*100;
plot(C_L,differ)


title('Difference of Drag as Lift Increases')
xlabel('Coefficent of Lift')
ylabel('Percent Difference')