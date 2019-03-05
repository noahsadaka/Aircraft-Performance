%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Author: Noah Sadaka
% Date: February 5 2019
% Course: Aerospace Vehicle Performance

% Purpose: Solve for Assignments 3, 4, and 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

% ASSIGNMENT NUMBER (3, 4, 5)
Q = 5;

% AIRCRAFT VARIABLES

% Aircraft Weight [lbs]
W=45000;
% Aircraft Wing Area [ft2]
S=520;
% ISA conditions [K]
ISA=10;
% Mean Aerodynamic Chord [ft]
MAC = 8.286;
% Tail Arm (ft}
LT = 40.56;
% Wing span [ft]
b = 67.85;
% Aspect Ratio
AR = b^2/S;
% Number of Engines
n_eng = 2;

% FLIGHT PROFILE

% Initial leg. VCAS1 in kt and ALT1 in ft
VCAS1 = 250;
ALT1 = 100*100;

% Second Leg. VCAS2 in kt
VCAS2 = 290;
ALT2 = 0;

%Third Leg. M3 in Mach, ALT4 is final altitude
M3 = 0.74;
ALT4 = 390*100;

% Datasheet Data Input

BUFF_M = [0.2750 0.3000 0.3250 0.3500 0.3750 0.4000 0.4250 0.4500 0.4750 0.5000 0.5250 0.5500 0.5750 0.6000 0.6250 0.6500 0.6750 0.7000 0.7250 0.7500 0.7750 0.8000 0.8250 0.8500 0.8750 0.9000];
BUFF_CL = [1.3424 1.3199 1.2974 1.2667 1.2310 1.1930 1.1551 1.1191 1.0863 1.0577 1.0337 1.0142 0.9989 0.9868 0.9764 0.9659 0.9530 0.9349 0.9085 0.8698 0.8149 0.7391 0.6373 0.5039 0.3330 0.118];

M_COMPR = [0.3 0.4 0.5 0.6 0.7 0.74 0.77 0.8 0.81 0.82 0.83 0.84 0.85];
CL_COMPR = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
CD_COMPR = [ 0 0 0 0 0 0 0 0;
    0 0 0 0 0.0001 0.0004 0.0013 0.003;
    0 0 0.0001 0.0004 0.0008 0.0017 0.0035 0.0072;
    0 0.0003 0.0006 0.001 0.0018 0.0033 0.0068 0.0134;
    0.0005 0.001 0.002 0.003 0.0042 0.0068 0.0114 0.0195;
    0.002 0.0025 0.0035 0.0045 0.0057 0.0084 0.0136 0.0225;
    0.003 0.0035 0.0045 0.0055 0.0066 0.0099 0.0157 0.026;
    0.004 0.0045 0.0055 0.0065 0.0081 0.0131 0.0214 0.0374;
    0.005 0.0055 0.0065 0.0075 0.00996 0.017 0.026 0.049;
    0.006 0.0065 0.0075 0.0095 0.013 0.021 0.031 0.065;
    0.007 0.0075 0.009 0.0115 0.017 0.026 0.036 0.08;
    0.0083 0.0095 0.0115 0.015 0.021 0.0315 0.043 0.1;
    0.01 0.0125 0.016 0.02 0.0271 0.0381 0.0554 0.1194 ];

% Load in Thrust Data
Thrust_Data_import = xlsread('5_ACPerf_data_thrust.xls');

% Create ISA+10 Data by averaging the ISA and ISA+20 values
% Also remove all the NaN values and make a single matrix

Thrust_data = zeros(round(size(Thrust_Data_import(:,1),1)/2,0),3);
Thrust_altitudes = zeros(round(size(Thrust_Data_import(:,1),1)/2,0),1);
Thrust_mach = zeros(round(size(Thrust_Data_import(:,1),1)/2,0),1);
index = 1;

for i=1:round(size(Thrust_Data_import(:,1),1)/2,0)
    for j=round(size(Thrust_Data_import(:,1),1)/2,0):round(size(Thrust_Data_import(:,1),1),0)
        if j~=i
            if Thrust_Data_import(i,1) == Thrust_Data_import(j,1) && Thrust_Data_import(i,2) == Thrust_Data_import(j,2)
                Thrust_altitudes(index) = Thrust_Data_import(i,1);
                Thrust_mach(index) = Thrust_Data_import(i,2);
                Thrust_data(index,1:3) = [Thrust_Data_import(i,1), Thrust_Data_import(i,2) , 0.5*(Thrust_Data_import(i,4) + Thrust_Data_import(j,4))];
                index=index+1;
            end
        end
    end
end

Thrust_altitudes = unique(Thrust_altitudes);
Thrust_mach = unique(Thrust_mach);
Thrust_mach = Thrust_mach(2:end);

Thrust = zeros(length(Thrust_mach),length(Thrust_altitudes));

for i=1:length(Thrust_mach)
    for j=1:length(Thrust_altitudes)
        for k=1:length(Thrust_data(:,1))
            if Thrust_data(k,1) == Thrust_altitudes(j) && Thrust_data(k,2) == Thrust_mach(i)
                Thrust(i,j) = Thrust_data(k,3);
            end
        end
    end
end

%%

%Aerodynamic Parameter Initialization

To = 288.15; % SL temperature [K]
Po = 2116.22; % SL pressure [lb/ft3]
ro=0.002377; % SL density [sl/ft3]
ao=661.47; % SL speed of sound [kt]

T = zeros(length(1:ALT4),1); % Temperature
P = zeros(length(1:ALT4),1); % Pressure
r = zeros(length(1:ALT4),1); % Density
theta = zeros(length(1:ALT4),1); % Temp Ratio
delta = zeros(length(1:ALT4),1); % Pressure Ratio
sigma = zeros(length(1:ALT4),1); % Density Ratio
a = zeros(length(1:ALT4),1); % Speed of Sound

% Assignment 3 Parameter Initialization

VC = zeros(length(1:ALT4),1);
VT = zeros(length(1:ALT4),1);
VE = zeros(length(1:ALT4),1);
M = zeros(length(1:ALT4),1);
CL = zeros(length(1:ALT4),1);
qc = zeros(length(1:ALT4),1);
hp = zeros(length(1:ALT4),1);
M_check = 0; % flag once transition altitude is reached

% Assignment 4 Parameter Initialization
Re = zeros(length(1:ALT4),1);
CDC = zeros(length(1:ALT4),1);
CD = zeros(length(1:ALT4),1);
CLCD = zeros(length(1:ALT4),1);
Nz = zeros(length(1:ALT4),1);

% Assignment 5 Parameter Initialization
D = zeros(length(1:ALT4),1);
RC = zeros(length(1:ALT4),1);
AF = zeros(length(1:ALT4),1); 
Th = zeros(length(1:ALT4),1); 

%%

% Calculations!

for i=1:ALT4
    h_p = i-1; % Pressure altitude
    hp(i) = h_p;
    
    % Atmospheric Parameter Calculations
    
    if h_p>36068 % if above tropopause
        T(i) = 0.7519*To+ISA;
        delta(i) = 0.22336*exp(-1*(h_p-36089)/20806);
        P(i) = delta(i).*Po;
        theta(i) = T(i)./To;
        sigma(i) = 0.29707*exp(-1*(h_p-36089)/20806);
        r(i) = sigma(i).*ro;
    else
        T(i) = (1 - h_p * 6.87535E-6)*To + ISA;
        delta(i) = (1 - h_p * 6.87535E-6)^5.2559;
        P(i) = delta(i).*Po;
        theta(i) = T(i)./To;
        sigma(i) = delta(i)/theta(i);
        r(i) = sigma(i).*ro;
    end
    a(i) = ao.*sqrt(theta(i));
    
    % Assignment 3 Flight Phase Calculations
    
    if h_p < ALT1
        VC(i) = VCAS1;
        qc(i) = Po*((1+0.2*(VC(i)./ao)^2)^3.5-1);
        M(i) = sqrt(5*((qc(i)./P(i)+1)^0.2857-1));
        VT(i) = a(i).*M(i);
        VE(i) = VT(i).*sqrt(sigma(i));
    else
        VC(i) = VCAS2;
        qc(i) = Po*((1+0.2*(VC(i)./ao)^2)^3.5-1);
        M(i) = sqrt(5*((qc(i)./P(i)+1)^0.2857-1));
        VT(i) = a(i).*M(i);
        VE(i) = VT(i).*sqrt(sigma(i));
        if M(i) > 0.74 || M_check == 1
            if M_check == 0 % Record the transition altitude and toggle the flag
                ALT2 = h_p;
                M_check = 1;
            end
            M(i) = 0.74;
            VT(i) = M(i).*a(i);
            VE(i) = VT(i).*sqrt(sigma(i));
            qc(i) = P(i).*(((M(i)^2/5)+1)^(1/0.2857)-1);
            VC(i) = ao*sqrt(5*((qc(i)./Po+1)^0.2857-1));
        end
    end
    CL(i) = (295.4*W)/(S*VE(i)^2);
    
    % Assignment 4 Calculations
    
    Re(i) = M(i) * 5.13384E6 * ((theta(i) + 0.38312)/(theta(i)^2))*delta(i); % Source: 11-10 in Boeing
    
    CDC(i) = interp2(CL_COMPR , M_COMPR , CD_COMPR , CL(i), M(i)); % Comp Drag
    
    CL_test = interp1(BUFF_M,BUFF_CL,M(i)); % CL at 9% MAC
    CL_FWD = CL_test * (1+ (MAC/LT) * (.25-.09)); % CL at 25% MAC
    Nz(i) = (1481.3*CL_FWD*M(i)^2*S*delta(i))/(W); % Buffet Load Facteur
    
    CD(i) = 0.0206 + .0364 * CL(i)^2+ CDC(i); % Drag Coeff, using values from perf data sheet and including compressibility
    
    CLCD(i) = CL(i)/CD(i); % lift to drag ratio
    
    % Assignment 5 Calculations
    
    Th(i) = n_eng * interp2(Thrust_altitudes,Thrust_mach,Thrust,hp(i),M(i));
    
    V = VT(i) * 1.6878; % true velocity in ft/s
    
    D(i) = 0.5 * CD(i) * r(i) * S * V^2;
    
    % Calculate AF
    if M_check == 0 % if below transition altitude
        % CONSTANT VCAS BELOW TROP
        fi = (((1 + 0.2*M(i)^2)^3.5) -1)/(0.7*M(i)^2 * (1+0.2*M(i)^2)^2.5);
        AF(i) = 0.7 * M(i)^2 * (fi - 0.190263 * ((T(i)-ISA)/T(i)));
        
    elseif M_check == 1 && hp(i) < 36089 % if between tran and trop
        % CONSTANT MACH BELOW TROP
        AF(i) = -0.133184*M(i)^2;
        
    else % if above tropopause
        % CONSTANT MACH ABOVE TROP
        AF(i) = 0;
    end
    
    
    
    
    RC(i) = V * 60 * ((Th(i)-D(i))/W)/(1+AF(i));

end

% Plot Solutions

if Q == 3
    
    figure(1)
    hold on
    yyaxis left
    plot(hp(:),VT(:))
    plot(hp(:),VC(:))
    yyaxis right
    plot(hp(:),M(:))
    plot(hp(:),CL(:))
    legend('VT','VC','M','CL')
    saveas(figure(1),'3_allplot.png')
    
elseif Q == 4
    
    figure(1)
    plot(hp(:), Re(:),'LineWidth',1.25)
    xlabel('Pressure Altitude [ft]','interpreter','latex','FontSize',13)
    ylabel('Reynolds Number','interpreter','latex','FontSize',13)
    title({'Reynolds Number at MAC v. Pressure Altitude' ' '},'interpreter','latex','FontSize',13)
    saveas(figure(1),'4_re.png')
    close(figure(1))
    
    figure(2)
    plot(hp(:), CDC(:),'LineWidth',1.25)
    xlabel('Pressure Altitude [ft]','interpreter','latex','FontSize',13)
    ylabel('Compressibility Drag Coefficient $C_{D,COMP}$','interpreter','latex','FontSize',13)
    title({'Compressibility Drag Coefficient v. Pressure Altitude' ' '},'interpreter','latex','FontSize',13)
    saveas(figure(2),'4_cdcomp.png')
    close(figure(2))
    
    figure(3)
    plot(hp(:), CD(:),'LineWidth',1.25)
    xlabel('Pressure Altitude [ft]','interpreter','latex','FontSize',13)
    ylabel('Total Drag Coefficient $C_D$','interpreter','latex','FontSize',13)
    title('Total Drag Coefficient v. Pressure Altitude','interpreter','latex','FontSize',13)
    saveas(figure(3),'4_cd.png')
    close(figure(3))
    
    figure(4)
    plot(hp(:), CLCD(:),'LineWidth',1.25)
    xlabel('Pressure Altitude [ft]','interpreter','latex','FontSize',13)
    ylabel('Lift to Drag Ratio $\frac{C_L}{C_D}$','interpreter','latex','FontSize',13)
    title('Lift to Drag Ratio v. Pressure Altitude','interpreter','latex','FontSize',13)
    saveas(figure(4),'4_clcd.png')
    close(figure(4))
    
    figure(5)
    plot(hp(:), Nz(:),'LineWidth',1.25)
    xlabel('Pressure Altitude [ft]','interpreter','latex','FontSize',13)
    ylabel('Buffet Load Factor $N_z$ [g]','interpreter','latex','FontSize',13)
    title('Load Factor at Buffet Onset v. Pressure Altitude','interpreter','latex','FontSize',13)
    saveas(figure(5),'4_nz.png')
    close(figure(5))
    
elseif Q == 5
    figure(1)
    plot(hp(:), RC(:),'LineWidth',1.25)
    xlabel('Pressure Altitude [ft]','interpreter','latex','FontSize',13)
    ylabel('Rate of Climb [ft/min]','interpreter','latex','FontSize',13)
    title('Rate of Climb v. Pressure Altitude','interpreter','latex','FontSize',13)
    saveas(figure(1),'5_rc.png')
    %close(figure(1))
    
    'at SL:'  %#ok<*NOPTS>
    RC(1)
    'at FL100 and 250KCAS'
    RC(10000)
    'at FL100 and 290 KCAS'
    RC(10001)
    'at Fl200'
    RC(20000)
    'at hp Tran'
    RC(ALT2)
    'after hp tran'
    RC(ALT2+1)
    'at fl330'
    RC(33000)
    ' at trop'
    RC(36089)
    'after trop'
    RC(36090)
    'at FL390'
    RC(end)
    
    'Question 2 answer'
    [~,Q2_ans] = min(abs(RC-300))
    Q2_ans
    
    % Calculate the time to climb 
    % SL to 10000
    
    RC_SL_10 = 0.5*(RC(1) + RC(9999));
    time_SL_10 = 10000/RC_SL_10;
    
    % Accelerate from 250 KCAS to 290 KCAS
    a_250 = (Th(9999)-D(9999))/W * 32.17;
    a_290 = (Th(10000)-D(10000))/W * 32.17;
    a_250_290 = 0.5*(a_250+a_290);
    time_10 = 1/60*(1.6878*(290-250)/a_250_290);
    
    % Climb from Fl 10 to FL 20
    RC_10_20 = 0.5*(RC(10000) + RC(10000));
    time_10_20 = 10000/RC_10_20;
    
    % Climb from FL 20 to Transition
    RC_20_tr = 0.5*(RC(20000) + RC(ALT2+1));
    time_20_tr = (ALT2-20000)/RC_20_tr;
    
    %Climb from tran to FL33
    RC_tr_33 = 0.5*(RC(ALT2+1) + RC(33000));
    time_tr_33 = (33000-ALT2)/RC_tr_33;
    
    total_climb_time = 60*(time_SL_10+time_10+time_10_20+time_20_tr+time_tr_33);
    
    ' Climb Time =' 
    total_climb_time
    
    
    
end
