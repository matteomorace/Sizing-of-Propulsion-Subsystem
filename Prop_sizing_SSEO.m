%% PROPULSIVE SYSTEM BREAKDOWN

clear all
clc

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% READ BEFORE: TO CHANGE FROM COMPUTATIONS WITH EXPECTED AND PREVIUSOLY
% COMPUTED VALUES OF DV IT IS NECESSARY TO COMMENT AND UNCOMMENT LINES
% 15/16, 26, 30-31/34-35, 50-59/62-70, 93/94, 151




%% Definition of quantities
DV = (1260.2); %expected case
% DV = 1.1*(1150.4); %required Delta V from computed manoeuvres from HW1
Is = 221.2;      %specific impulse (average of Is during the manoeuvres) REFERENCE FIL
DM = 949;
DM_mrgin = 1.2*949;   %dry mass with 20% system margin   FROM REFERENCE!!!!!!!!!!!!!!!!

g = 9.81;

MR = exp(DV/(Is*g)); % mass ratio m0/mf Tsiolkovsky
M_tot = MR*DM_mrgin;       %total mass at launch
Mprop = (M_tot-DM_mrgin);
% Mprop = 1.055*Mprop;   % to uncomment only in case of computed case sizing
M_tot = Mprop + DM_mrgin;
% Prop mass calculations at each manoeuvre

Is = [219.7 225 225 225 225 225 225 224 224 224]; %expected case
DV = [28 567 98 102.8 113 37.4 60 150 65 39]; %expected case


% Is = [219.7 225 225 225 225 225 224 224 224]; %computed case
% DV = [28 524.1 214.6 91 40.1 41.1 146.5 65 40]; %computed case


Mi = M_tot;
M_prop_consumed = zeros(1,length(Is));
for ii = 1:length(Is)
    I = Is(ii);
    Dv = DV(ii);
    MR = exp(Dv/(I*g));
    Mf = Mi/MR;
    mprop_consumed = Mi-Mf;
    Mi = Mf;
    M_prop_consumed(ii) = mprop_consumed;
end
% computed with data from reference (expected DV)
fprintf("Propellant mass consumed during MCC is: %f kg \n", M_prop_consumed(1));
fprintf("Propellant mass consumed during LOI1 is: %f kg \n", M_prop_consumed(2));
fprintf("Propellant mass consumed during LOI2 is: %f kg \n", M_prop_consumed(3));
fprintf("Propellant mass consumed during LOI3 is: %f kg \n", M_prop_consumed(4));
fprintf("Propellant mass consumed during LOI4 is: %f kg \n", M_prop_consumed(5));
fprintf("Propellant mass consumed during LOI5 is: %f kg \n", M_prop_consumed(6));
fprintf("Propellant mass consumed during MOI is: %f kg \n", M_prop_consumed(7));
fprintf("Propellant mass consumed during SK is: %f kg \n", M_prop_consumed(8));
fprintf("Propellant mass consumed for Margin is: %f kg \n", M_prop_consumed(9));
fprintf("Propellant mass consumed during Extended mission is: %f kg \n", M_prop_consumed(10));

% computed with data our computations (computed DV)
% fprintf("Propellant mass consumed during MCC is: %f kg \n", M_prop_consumed(1));
% fprintf("Propellant mass consumed during LOI1 is: %f kg \n", M_prop_consumed(2));
% fprintf("Propellant mass consumed during LOI3 is: %f kg \n", M_prop_consumed(3));
% fprintf("Propellant mass consumed during LOI4 is: %f kg \n", M_prop_consumed(4));
% fprintf("Propellant mass consumed during LOI5 is: %f kg \n", M_prop_consumed(5));
% fprintf("Propellant mass consumed during MOI is: %f kg \n", M_prop_consumed(6));
% fprintf("Propellant mass consumed during SK is: %f kg \n", M_prop_consumed(7));
% fprintf("Propellant mass consumed for Margin is: %f kg \n", M_prop_consumed(8));
% fprintf("Propellant mass consumed during Extended mission is: %f kg \n", M_prop_consumed(9));


%% Tank sizing

rho = 1008; %density of hydrazine

Vprop = (Mprop/rho); % propellant volume
V_tank = 1.2*Vprop; % tank volume with 20% margin



P_tank_pressurefed = 18.6; % tank pressure when pressure fed activated FROM REFERENCE arc3!!!!!!!!!!!!!!!!  [bar]
r_tank = ((3*V_tank/2)/(4*pi))^(1/3); % radius of one tank considering 2 propellant tanks
sigma = 880*10^6;  % tensile strength of titanium alloy Ti6A14V
t_tank = P_tank_pressurefed*10^5*r_tank/(2*sigma);

m_tank = 4400*(4/3)*pi*((r_tank+t_tank)^3-r_tank^3) 

%% BLOWDOWN PART

V_blowdown_in = V_tank - Vprop; % initial blowdown gas volume

Mprop_MCC = 26.4; %propellant mass expected used during MCC manoeuvre FROM REFERENCE!!!!!!!!!!!!!!!!
% Mprop_MCC = 26.97; %computed case
Vprop_post_MCC = Vprop - (Mprop_MCC/rho); %propellant volume after MCC
V_blowdown_out = V_tank - Vprop_post_MCC; %final blowdown gas volume

B = V_blowdown_out/V_blowdown_in; %blowdown ratio

Ptank = 14.1;  %propellant tank pressure FROM REFERENCE arc3!!!!!!!!!!!!!!!!  [bar]
P_blowdown_out = Ptank;  %final blowdown gas pressure
P_blowdown_in = P_blowdown_out*B; %initial blowdown gas pressure [bar]

mol_blowdown = (P_blowdown_in*10^5*V_blowdown_in)/(8.31446261815324*298);
m_blowdown = mol_blowdown*4.003e-3;


%% VARIABLES FOR THE PRESSURIZER

T_tank = 298;                                       % Temperature of all gas and fuel
gamma = 1.67;                                       % constant of helium
R = 2077.3;                                         % [J/kgK]
M_mol = 4.003e-3;                                   % Molar mass  kg/mol

P_pressurizer = P_tank_pressurefed*15*10^5;                                              % Pressurizer pressure 
mol_pressureant_nominal = (P_tank_pressurefed*10^5*V_tank)/(8.31446261815324*T_tank)-mol_blowdown;    % Numero di moli of the pressurizer
m_pressurizer = M_mol*mol_pressureant_nominal;                                           % Mass of pressurizer 
m_pressurizer_margin = m_pressurizer * 1.2;                                              % Mass of pressurizer with 20% margin
V_pres = mol_pressureant_nominal*8.314462618*T_tank/P_pressurizer;                       % Volume pressurizer
V_pres_margin = 1.1*V_pres;     % Tank Volume with 10% margin
h_pressuriser = 0.752; %from reference
r_pressurizer = sqrt(V_pres/(pi*h_pressuriser));
r_pressurizer_margin = sqrt(V_pres_margin/(pi*h_pressuriser));                           % radius of one tank considering 2 propellant tanks
sigma = 950e6;                                                                           % tensile strength of titanium alloy Ti6A14V
t_pressurizer = P_pressurizer*r_pressurizer/(sigma);                                     % Thickness of pressurizer tank

m_tank_p = 2780*(pi*(r_pressurizer+t_pressurizer)^2*(h_pressuriser+2*t_pressurizer)-pi*r_pressurizer^2*h_pressuriser)  %VEDERE ESERCITAZIONE


%% Pressure at each manoeuvre


V_tank_initial = V_tank;
Mprop_initial = Mprop;
P_sequence = zeros(1,length(M_prop_consumed)-3);

for ii = 1:length(M_prop_consumed)-3
    V_tank_left = V_tank_initial - (Mprop_initial-M_prop_consumed(ii+1))/rho;
    mol_He = (P_tank_pressurefed*10^5*V_tank_left)/(8.31446261815324*T_tank)-mol_blowdown;
    P_pressurizer_man = ((mol_pressureant_nominal-mol_He)*8.314462618*T_tank/V_pres_margin)/10^5;
    P_sequence(ii) = P_pressurizer_man;
    Mprop_initial = Mprop_initial - M_prop_consumed(ii+1);
end
% expected case
fprintf("Pressure in press. tank after LOI1 is: %f bar \n", P_sequence(1));
fprintf("Pressure in press. tank after LOI2 is: %f bar \n", P_sequence(2));
fprintf("Pressure in press. tank after LOI3 is: %f bar \n", P_sequence(3));
fprintf("Pressure in press. tank after LOI4 is: %f bar \n", P_sequence(4));
fprintf("Pressure in press. tank after LOI5 is: %f bar \n", P_sequence(5));
fprintf("Pressure in press. tank after MOI is: %f bar \n", P_sequence(6));
fprintf("Pressure in press. tank after SK is: %f bar \n", P_sequence(7));
% comment line 151 if computed case


%% Power and Mass budget
M_budget = 1.1*(2*m_tank + m_tank_p + m_pressurizer_margin + 8*0.59 + 4*1.12) ;
M_budget_actual = 1.1*(2*35.4+ 12.7 +3.3 + 8*0.59 + 4*1.12) ;
