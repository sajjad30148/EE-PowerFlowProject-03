%% Assignment 03 - Continuation Power Flow - Sajjad Uddin Mahmud - Fall 2022 - WSU

%% Basic Initialization
clc;
clear all;
close all;

%% Reading From Bus Data
%% Bus Number
All_Bus_Number = xlsread('IEEE14_Formatted.xlsx','A3:A16'); % Reading All Bus ID Data
Total_Bus = length(All_Bus_Number); % Calculating Total Bus Number

%% Bus Type
All_Bus_Type = xlsread('IEEE14_Formatted.xlsx','G3:G16'); % Reading All Bus Type Data
PQ_Bus_Type = 0; 
PQ_Bus_Type_1 = 1; 
PV_Bus_Type = 2; 
Slack_Bus_Type = 3;
Bus = All_Bus_Number(find(All_Bus_Type ~= Slack_Bus_Type)); % Bus Type Data Except the Slack Bus
Bus_Type = All_Bus_Type(Bus); % Bus Type Data Except the Slack Bus

%% Bus Information
% Slack_Bus_Number = 1

Base_MVA = 100; 
V_Desired = xlsread('IEEE14_Formatted.xlsx','O3:O16'); % Given Initial Voltage
P_Load_All = xlsread('IEEE14_Formatted.xlsx','J3:J16')/Base_MVA; % Load MW pu
P_Load = P_Load_All(find(All_Bus_Type~=Slack_Bus_Type ));
Q_Load_All = xlsread('IEEE14_Formatted.xlsx','K3:K16')/Base_MVA; % Load MVAR pu
Q_Load = Q_Load_All(find(All_Bus_Type~=Slack_Bus_Type ));
P_Gen_All = xlsread('IEEE14_Formatted.xlsx','L3:L16')/Base_MVA; % Generator MW pu
P_Gen = P_Gen_All(find(All_Bus_Type~=Slack_Bus_Type ));
Q_Gen_All = xlsread('IEEE14_Formatted.xlsx','M3:M16')/Base_MVA; % Generator MVAR pu
Q_Gen = Q_Gen_All(find(All_Bus_Type~=Slack_Bus_Type ));

%% Reading from Branch Data
%% Branch Number
To_Bus = xlsread('IEEE14_Formatted.xlsx','A19:A38'); 
From_Bus = xlsread('IEEE14_Formatted.xlsx','B19:B38');
Total_Bus = max(max(To_Bus),max(From_Bus));

%% Bus Shunt Conductance and Shunt Susceptance
G_Shunt_Bus = xlsread('IEEE14_Formatted.xlsx','R3:R16');
B_Shunt_Bus = xlsread('IEEE14_Formatted.xlsx','S3:S16');

%% Calculating Bus Shunt Admittance
Y_Shunt_Bus = G_Shunt_Bus + j.*B_Shunt_Bus;

%% Branch Resistance Per Unit
R_Branch = xlsread('IEEE14_Formatted.xlsx','G19:G38');

%% Branch Reactance Per Unit
X_Branch = xlsread('IEEE14_Formatted.xlsx','H19:H38');

%% Line Charging B Per Unit
B_Branch = xlsread('IEEE14_Formatted.xlsx','I19:I38');

%% Transformer Turns Ratio
XFR_TurnRatio = xlsread('IEEE14_Formatted.xlsx','O19:O38');

%% Calculating Branch Impedence and Admittance
for i=1:length(From_Bus)
    Z_Branch(i) = R_Branch(i) + j*X_Branch(i); % Per Unit Impedance
    Y_Branch(i) = 1/Z_Branch(i); % Per Unit Admittance
end 

%% Tap Consideration
Tap_Consideration = 1; % 0 = Without Taps, 1 = With Taps
if (Tap_Consideration == 0)
    for i = 1:length(XFR_TurnRatio)
        XFR_TurnRatio(i) = 0; % If We Do Not Consider Tap, All the Turn Ratio of Transformer are 0
    end
end

%% Calculating Y Bus Matrix: 

% Initialization
Y_Bus = zeros(Total_Bus,Total_Bus);

% LOOP: Computing Off-Diagonal Elements
for i=1:length(Y_Branch)
    if (XFR_TurnRatio(i)==0)
        Y_Bus(To_Bus(i),From_Bus(i)) = - Y_Branch(i);
        Y_Bus(From_Bus(i),To_Bus(i)) = - Y_Branch(i);
        Y_Bus_Diag(To_Bus(i),From_Bus(i)) = - Y_Branch(i);
        Y_Bus_Diag(From_Bus(i),To_Bus(i)) = - Y_Branch(i);
    else
        T = (1/(XFR_TurnRatio(i)));
        Y_Bus(To_Bus(i),From_Bus(i)) = - Y_Branch(i) * (T);
        Y_Bus(From_Bus(i),To_Bus(i)) = - Y_Branch(i) * (T);
        Y_Bus_Diag(To_Bus(i),From_Bus(i)) = - Y_Branch(i);
        Y_Bus_Diag(From_Bus(i),To_Bus(i)) = - Y_Branch(i) * (T^2);
    end
end

% LOOP: Computing Diagonal Elements
Y_Bus_Sum = sum(Y_Bus_Diag);
for i=1:Total_Bus
    Y_Bus(i,i) = - Y_Bus_Sum(i) + Y_Shunt_Bus(i); % Adding Shunt Capacitance
end

% LOOP: Adding Line Charaging Capacitance
for i=1:length(To_Bus)
    Y_Bus(To_Bus(i),To_Bus(i)) = Y_Bus(To_Bus(i),To_Bus(i)) + j * (B_Branch(i) / 2);
    Y_Bus(From_Bus(i),From_Bus(i)) = Y_Bus(From_Bus(i),From_Bus(i)) + j * (B_Branch(i) / 2);
end

% Converting Y Bus Data into Polar Form
Rho = abs(Y_Bus); % Magnitude of Y Bus Entries
Theta = angle(Y_Bus); % Angle of Y Bus Entries in radian
B = imag(Y_Bus); % Imaginary Part of Y Bus Entries
G = real(Y_Bus); % Real Part of Y Bus Entries

% End of Y Bus Formation. Y Bus is Ready


%% Getting K Vector
PQ_Bus_Number = All_Bus_Number(find((All_Bus_Type ~= Slack_Bus_Type) & (All_Bus_Type ~= PV_Bus_Type)));
PQ_Bus_Title = All_Bus_Number(find((All_Bus_Type ~= Slack_Bus_Type) & (All_Bus_Type ~= PV_Bus_Type)));

P_Scheduled = P_Gen_All - P_Load_All;
Q_Scheduled = Q_Gen_All-Q_Load_All;
P_K = P_Scheduled(find(All_Bus_Type ~= Slack_Bus_Type));
Q_K = Q_Scheduled(find((All_Bus_Type~=Slack_Bus_Type) & (All_Bus_Type~=PV_Bus_Type)));
K_Vector = [P_K ; Q_K];


for Parameter = 1:length(PQ_Bus_Title)

    %% Lambda and Sigma Initialization
    Lambda = 0;
    Sigma = 0.1;
    Bus_Nose = PQ_Bus_Title(Parameter);
    
    % Initial Voltage Magnitude and Angle (Flat Start)
    
    V = ones(1,length(All_Bus_Number));
    Delta_in_Rad = zeros(1,length(All_Bus_Number));
    
    %% Newton-Raphson Power Flow Solution
    Tolerance = 0.01;
    
    % Function Calling: Power Flow Using Newton-Raphson Method
    [V, Delta_in_Rad, Iteration] = Newton_Raphson_Function(Y_Bus, V, Delta_in_Rad, Lambda, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type, Tolerance);
    
    Nose = 0;
    
    while 1
    
        Nose = Nose + 1;
        V_X = V;
        Delta_X = Delta_in_Rad;
        Lambda_X = Lambda;
        
        %% Upper Part
        
        % Predictor
        [J] = Jacobian_Function(V, Delta_in_Rad, Y_Bus, All_Bus_Type);
        Delta_Base = Delta_in_Rad(find(All_Bus_Type ~= Slack_Bus_Type ));
        V_Base = V(find((All_Bus_Type ~= Slack_Bus_Type) & (All_Bus_Type ~= PV_Bus_Type)));
        
        % Initializing Extra Equation for Lambda
        E_K_Vector = zeros(1, length(K_Vector)+1);
        E_K_Vector(end) = 1; % Setting Last Value of Extra Vector as 1
        
        % Initializing Predictor B Vector
        Predictor_B_Vector = zeros(length(E_K_Vector),1);
        Predictor_B_Vector(end) = 1;
        
        % Calculating Tangent Vector Using LU Factorization Dolittle's Method
        Delta_V_Lamda_Vector = [transpose(Delta_Base); transpose(V_Base); Lambda];
        J_Predictor = [J K_Vector; E_K_Vector];
        [Tangent_Vector] = LU_Factorization_Dolittle_Function(J_Predictor, Predictor_B_Vector);
        
        % Computing Predicted Values
        Delta_V_Lamda_Vector_Predicted = Delta_V_Lamda_Vector + (Sigma * Tangent_Vector);
        
        % Predicted Lambda
        Lambda = Delta_V_Lamda_Vector_Predicted(end);
        
        % LOOP: Predicted Voltages and Angles
        for i=1:(length(Delta_V_Lamda_Vector_Predicted) - 1)
            if (i<=length(Delta_Base))
                Delta_Predicted(i) = Delta_V_Lamda_Vector_Predicted(i);
            else
                V_Predicted(i-(length(Delta_Base))) = Delta_V_Lamda_Vector_Predicted(i);
            end
        end
        
        % Preparing for Correction Stage
        V_Index = (find((All_Bus_Type == PQ_Bus_Type) | (All_Bus_Type == PQ_Bus_Type_1)));
        Delta_Index = find(All_Bus_Type ~= Slack_Bus_Type );
        
        for i=1:length(Delta_Index)
            Delta_in_Rad(Delta_Index(i)) = Delta_Predicted(i);
        end
        
        for i=1:length(V_Index)
            V(V_Index(i)) = V_Predicted(i);
        end
        
        % Corrector
        
        % Function Calling: Power Flow Using Newton-Raphson Method
        [V, Delta_in_Rad, Iteration] = Newton_Raphson_Function(Y_Bus, V, Delta_in_Rad, Lambda, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type, Tolerance);
        
        if Iteration >= 25
            V = V_X;
            Delta_in_Rad = Delta_X;
            Lambda = Lambda_X;
            break;
        end
        
        V_PV_Upper(Nose) = V(Bus_Nose);
        Delta_PV_Upper(Nose) = Delta_in_Rad(Bus_Nose);
        Lamda_PV_Upper(Nose) = Lambda;
    
    end
    
    %% Critial Point

    Sigma = 0.025; % Changing Sigma to 0.025 Because Voltage Changes at a Slower Rate than Power
    Critical = 0;
    
    while 1  
    
        Critical = Critical + 1;
            
        V_X = V;
        Delta_X = Delta_in_Rad;
        Lambda_X = Lambda;
        
        % Predictor Stage
        
        % Function Calling: Jacobian
        [J] = Jacobian(V, Delta_in_Rad, Y_Bus, All_Bus_Type);
        
        Delta_Base = Delta_in_Rad(find(All_Bus_Type ~= Slack_Bus_Type ));
        V_Base = V(find((All_Bus_Type ~= Slack_Bus_Type) & (All_Bus_Type ~= PV_Bus_Type)));
        
        % Initializing Extra Equation for Lambda
        E_K_Vector = zeros(1, length(K_Vector) + 1);
        E_K_Vector(Total_Bus - 1 + find(PQ_Bus_Number == Bus_Nose)) = -1;
        
        % Initializing Predictor B Vector
        Predictor_B_Vector = zeros(length(E_K_Vector),1);
        Predictor_B_Vector(end) = 1;
        
        % Calculating Tangent Vector Using LU Factorization Dolittle's Method
        Delta_V_Lamda_Vector = [transpose(Delta_Base); transpose(V_Base); Lambda];
        J_Predictor = [-J K_Vector; E_K_Vector];
        [Tangent_Vector] = LU_Factorization_Dolittle_Function(J_Predictor, Predictor_B_Vector);
        
        % Computing Predicted Values
        Delta_V_Lamda_Vector_Predicted = Delta_V_Lamda_Vector + (Sigma * Tangent_Vector);
        
        % Predicted Lambda
        Lambda = Delta_V_Lamda_Vector_Predicted(end);
        
        % LOOP: Predicted Voltages and Angles
        for i=1:(length(Delta_V_Lamda_Vector_Predicted) - 1)
            if (i <= length(Delta_Base))
                Delta_Predicted(i) = Delta_V_Lamda_Vector_Predicted(i);
            else
                V_Predicted(i-(length(Delta_Base))) = Delta_V_Lamda_Vector_Predicted(i);
            end
        end
        
        % Preparing for Correction Stage
        V_Index = (find((All_Bus_Type == PQ_Bus_Type) | (All_Bus_Type == PQ_Bus_Type_1)));
        Delta_Index = find(All_Bus_Type ~= Slack_Bus_Type );
        
        for i=1:length(Delta_Index)
            Delta_in_Rad(Delta_Index(i)) = Delta_Predicted(i);
        end
        
        for i=1:length(V_Index)
            V(V_Index(i)) = V_Predicted(i);
        end
        
        % Corrector Stage
        
        % Initialization
        Iteration = 0; 
        while 1
        
            Iteration = Iteration + 1;    
            PV = find(Bus_Type == PV_Bus_Type);
            
            % Schedule Real and Reactive Power
            P_Scheduled = transpose(P_Gen - (Lambda * P_Load));
            Q_Scheduled = transpose(Q_Gen - (Lambda * Q_Load));
            
            %% Calculating Real Power 
            
            % Initialization
            P_Calculated = zeros(1,length(Bus));
            
            % LOOP: Computing Real Power
            for i=1:length(Bus)
                for n=1:Total_Bus
                    P_Calculated(i) = P_Calculated(i) + (abs(abs(Y_Bus(Bus(i),n)) * V(Bus(i)) * V(n))) * (cos(angle(Y_Bus(Bus(i),n)) + Delta_in_Rad(n) - Delta_in_Rad(Bus(i))));
                end
            end
            
            %% Calculating Reactive Power 
            
            % Initialization
            Q_Calculated=zeros(1,length(Bus));
            
            % LOOP: Computing Reactive Power    
            for i=1:length(Bus)
                for n=1:Total_Bus
                    Q_Calculated(i) = Q_Calculated(i) + (abs(abs(Y_Bus(Bus(i),n)) * V(Bus(i)) * V(n))) * (sin(angle(Y_Bus(Bus(i),n)) + Delta_in_Rad(n) - Delta_in_Rad(Bus(i))));
                end
                Q_Calculated(i) = - Q_Calculated(i);
            end
            
            %% Calculating Mismatch
            Delta_P = P_Scheduled - P_Calculated;
            Delta_Q = Q_Scheduled - Q_Calculated;
            Delta_Q(:,PV)=[];
            Delta_P_Q_0 = [transpose(Delta_P);transpose(Delta_Q); 0];
            
            % Initializing Extra Equation for Lambda
            E_K_Vector = zeros(1, length(K_Vector)+1);
            E_K_Vector(Total_Bus - 1 + find(PQ_Bus_Number == Bus_Nose)) = -1;
            
            % Function Calling: Jacobian
            [J] = Jacobian_Function(V,Delta_in_Rad, Y_Bus, All_Bus_Type);
                            
            J_Corrected = [-J K_Vector; E_K_Vector];

            Mismatch_V_Delta_Lambda = LU_Factorization_Dolittle_Function(J_Corrected, -Delta_P_Q_0);
            Delta_Lamda = Mismatch_V_Delta_Lambda(end);
            Lambda = Lambda + Delta_Lamda;
            
            for i=1:(length(Mismatch_V_Delta_Lambda) - 1)
                if (i <= length(Delta_Index))
                    Mismatch_Delta(i) = Mismatch_V_Delta_Lambda(i);
                else
                    Mismatch_V(i-(length(Delta_Index))) = Mismatch_V_Delta_Lambda(i);
                end
            end
            
            for i=1:length(Delta_Index)
                Delta_in_Rad(Delta_Index(i)) = Mismatch_Delta(i) + Delta_in_Rad(Delta_Index(i));
            end
            
            for i=1:length(V_Index)
                V(V_Index(i)) = Mismatch_V(i) + V(V_Index(i));
            end
            
            if (max(abs(Mismatch_V_Delta_Lambda)) < 0.1)
                break;
            end
            
            if Iteration >= 5
                break;
            end
        end
        
        if Iteration >= 5
            V = V_X;
            Delta_in_Rad = Delta_X;
            Lambda = Lambda_X;
            break;
        end
        
        V_PV_Critical(Critical) = V(Bus_Nose);
        Delta_PV_Critical(Critical) = Delta_in_Rad(Bus_Nose);
        Lamda_PV_Critical(Critical) = Lambda;
    end
    
    %% Lower Part
    Sigma = 0.1; 
    
    Lower = 0;

    while 1
        Lower = Lower + 1;
        V_X = V;
        Delta_X = Delta_in_Rad;
        Lambda_X = Lambda;
        
        % Predictor Stage
        [J] = Jacobian(V,Delta_in_Rad, Y_Bus, All_Bus_Type);
        Delta_Base = Delta_in_Rad(find(All_Bus_Type ~= Slack_Bus_Type ));
        V_Base = V(find((All_Bus_Type ~= Slack_Bus_Type) & (All_Bus_Type ~= PV_Bus_Type)));
        
        % Initializing Extra Equation for Lambda
        E_K_Vector = zeros(1, length(K_Vector) + 1);
        E_K_Vector(length(E_K_Vector)) = -1;

        % Initializing Predictor B Vector
        Predictor_B_Vector = zeros(length(E_K_Vector),1);
        Predictor_B_Vector(length(E_K_Vector)) = 1;
        
        % Calculating Tangent Vector Using LU Factorization Dolittle's Method
        Delta_V_Lamda_Vector = [transpose(Delta_Base); transpose(V_Base); Lambda];
        J_Predictor = [J K_Vector; E_K_Vector];
        [Tangent_Vector] = LU_Factorization_Dolittle_Function(J_Predictor, Predictor_B_Vector);
        Delta_V_Lamda_Vector_Predicted = Delta_V_Lamda_Vector + (Sigma * Tangent_Vector);

        % Predicted Lambda
        Lambda = Delta_V_Lamda_Vector_Predicted(end);
        
        % LOOP: Predicted Voltages and Angles
        for i=1:(length(Delta_V_Lamda_Vector_Predicted)-1)
            if (i<=length(Delta_Base))
                Delta_Predicted(i)=Delta_V_Lamda_Vector_Predicted(i);
            else
                V_Predicted(i-(length(Delta_Base)))=Delta_V_Lamda_Vector_Predicted(i);
            end
        end
        
        % Preparing for Correction Stage
        V_Index=(find((All_Bus_Type == PQ_Bus_Type) | (All_Bus_Type == PQ_Bus_Type_1)));
        Delta_Index=find(All_Bus_Type ~= Slack_Bus_Type );
        
        for i=1:length(Delta_Index)
            Delta_in_Rad(Delta_Index(i)) = Delta_Predicted(i);
        end
        
        for i=1:length(V_Index)
            V(V_Index(i))=V_Predicted(i);
        end
        
        % Corrector Stage
        [V, Delta_in_Rad, Iteration] = Newton_Raphson_Function(Y_Bus, V, Delta_in_Rad, Lambda, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type, Tolerance);
        
        if Lambda <= 0.1
            V = V_X;
            Delta_in_Rad = Delta_X;
            Lambda = Lambda_X;
            break;
        end
        
        V_PV_Lower(Lower) = V(Bus_Nose);
        Delta_PV_Lower(Lower) = Delta_in_Rad(Bus_Nose);
        Lamda_PV_Lower(Lower) = Lambda;
    end
    
    V_PV = [V_PV_Upper V_PV_Critical V_PV_Lower];
    Lamda_PV = [Lamda_PV_Upper Lamda_PV_Critical Lamda_PV_Lower];
    

    %% Output

    figure (Parameter)
    plot (Lamda_PV, V_PV,'--o')
    hold on;
    plot(Lamda_PV_Upper, V_PV_Upper,'g')
    hold on;
    plot(Lamda_PV_Critical, V_PV_Critical,'r')
    hold on
    plot(Lamda_PV_Lower, V_PV_Lower,'b')
    
    
    title(['PV Curve of Bus ', num2str(Bus_Nose)])
    xlabel('Lamda');
    ylabel('Voltage Magnitude');
    
    V_PV = [];
    Lamda_PV = [];
    V_PV_Upper = [];
    V_PV_Critical = [];
    V_PV_Lower = [];
    Lamda_PV_Upper = [];
    Lamda_PV_Critical = [];
    Lamda_PV_Lower = [];
    Lambda = 0;

end