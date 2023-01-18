%% Power Flow Function: Newton-Raphson 

function [V, Delta_in_Rad, Iteration] = Newton_Raphson_Function(Y_Bus, V, Delta_in_Rad, Lambda, P_Gen, P_Load, Q_Gen, Q_Load, All_Bus_Number, All_Bus_Type, Tolerance)

%% Bus Number
Total_Bus = length(All_Bus_Number);

%% Bus Type
PQ_Bus_Type = 0; 
PQ_Bus_Type_1 = 1; 
PV_Bus_Type = 2; 
Slack_Bus_Type = 3;
Bus = All_Bus_Number(find(All_Bus_Type ~= Slack_Bus_Type)); % Bus Type Data Except the Slack Bus
Bus_Type = All_Bus_Type(Bus); % Bus Type Data Except the Slack Bus
PV_Bus = find(Bus_Type == PV_Bus_Type);

%% Initialization
Iteration = 0;
Tolerance = Tolerance;
while 1

Iteration = Iteration+1;  

%% Schedule Real and Reactive Power
P_Scheduled = transpose(P_Gen - (Lambda * P_Load));
Q_Scheduled = transpose(Q_Gen - (Lambda * Q_Load));

%% Calculating Real Power

% Initialization
P_Calculated = zeros(1,length(Bus));

% LOOP: Computing Real Power
for i=1:length(Bus)
    for n=1:Total_Bus
        P_Calculated(i) = P_Calculated(i) + (abs(abs(Y_Bus(Bus(i),n)) * V(Bus(i))*V(n))) * (cos(angle(Y_Bus(Bus(i),n)) + Delta_in_Rad(n) - Delta_in_Rad(Bus(i))));
    end
end

%% Calculating Reactive Power

% Initialization
Q_Calculated = zeros(1,length(Bus));


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
Delta_Q(:,PV_Bus)=[];
Delta_P_Q = [transpose(Delta_P);transpose(Delta_Q)];

[J] = Jacobian_Function(V, Delta_in_Rad, Y_Bus, All_Bus_Type);
Delta_J = Delta_in_Rad(find(All_Bus_Type ~= Slack_Bus_Type ));
V_J = V(find((All_Bus_Type ~= Slack_Bus_Type) & (All_Bus_Type ~= PV_Bus_Type)));


%% Updating V and Delta through LU Factorization

% Function Calling: LU Factorization Using Dolittle's Method
[V_Delta_Corrected] = LU_Factorization_Dolittle_Function(J,Delta_P_Q);
V_Delta_Corrected = transpose(V_Delta_Corrected);

% LOOP: Sorting the Voltages and Angles after LU Factorization
for i=1:length(V_Delta_Corrected)
    if (i <= length(Delta_P))
        Delta_Corrected(i) = V_Delta_Corrected(i);
    else
        V_Corrected(i-length(Delta_P)) = V_Delta_Corrected(i);
    end
end

% Updating Voltages and Angles
Delta_Updated = Delta_J + Delta_Corrected;
V_Updated = V_J .* (1 + V_Corrected);

% Preparing for Next Iteration
V_i = (find((All_Bus_Type == PQ_Bus_Type) | (All_Bus_Type == PQ_Bus_Type_1)));
Delta_i = find(All_Bus_Type ~= Slack_Bus_Type);

for i=1:length(Delta_i)
    Delta_New(Delta_i(i)) = Delta_Updated(i);
end

for i=1:length(V_i)
    V(V_i(i)) = V_Updated(i);
end

Delta_in_Rad = Delta_New;
Delta_in_Degree = (180 / pi) * Delta_in_Rad;

if (max(abs(Delta_P)) < Tolerance & max(abs(Delta_Q)) < Tolerance)
    break;
end

if Iteration >= 25
    break;
    
end
end
end

