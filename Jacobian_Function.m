%% Jacobian Function 

function [J] = Jacobian_Function(V, Delta_in_Rad, Y_Bus, All_Bus_Type)

PQNum=0;
PVNum=0;
countPQ=0;
countPV=0;
nbus=length(All_Bus_Type);
for i=1:nbus
    if All_Bus_Type(i) == 0
        countPQ=countPQ+1;
        PQNum=PQNum+1;
        PQBus(countPQ)=i;
    elseif All_Bus_Type(i) == 1
        countPQ=countPQ+1;
        PQNum=PQNum+1;
        PQBus(countPQ)=i;
    elseif All_Bus_Type(i) == 2
        countPV=countPV+1;
        PVNum=PVNum+1;
        PVBus(countPV)=i;
    end
end

% Calculating J11
J11 = zeros(nbus-1,nbus-1);
for i=1:nbus-1
    n=i+1;
    for jj=1:nbus-1
        m=jj+1;
        if n==m
            for k=1:nbus
                J11(i,jj) = J11(i,jj) + V(n)*V(k)*abs(Y_Bus(n,k))*sin(angle(Y_Bus(n,k))+Delta_in_Rad(k)-Delta_in_Rad(n));
            end
            J11(i,jj) = J11(i,jj) - V(n)^2*imag(Y_Bus(n,n));
        else
            J11(i,jj) = -V(n)*V(m)*abs(Y_Bus(n,m))*sin(angle(Y_Bus(n,m))+Delta_in_Rad(m)-Delta_in_Rad(n));
        end
    end
end

% Calculating J21
J21 = zeros(nbus-PVNum-1,nbus-1);
for i=1:nbus-PVNum-1
    n=PQBus(i);
    for jj=1:nbus-1
        m=jj+1;
        if n==m
            for k=1:nbus
                J21(i,jj) = J21(i,jj) + V(n)*V(k)*abs(Y_Bus(n,k))*cos(angle(Y_Bus(n,k))+Delta_in_Rad(k)-Delta_in_Rad(n));
            end
            J21(i,jj) = J21(i,jj) - V(n)^2*real(Y_Bus(n,n));
        else
            J21(i,jj) = -V(n)*V(m)*abs(Y_Bus(n,m))*cos(angle(Y_Bus(n,m))+Delta_in_Rad(m)-Delta_in_Rad(n));
        end
    end
end
% Calculating J12
J12 = zeros(nbus-1,nbus-PVNum-1);
for i=1:nbus-1
    n=i+1;
    for jj=1:nbus-PVNum-1
        m=PQBus(jj);
        if n==m
            for k=1:nbus
                J12(i,jj) = J12(i,jj) + V(n)*V(k)*abs(Y_Bus(n,k))*cos(angle(Y_Bus(n,k))+Delta_in_Rad(k)-Delta_in_Rad(n));
            end
            J12(i,jj) = J12(i,jj) + V(n)^2*real(Y_Bus(n,n));
        else
            J12(i,jj) = V(n)*V(m)*abs(Y_Bus(n,m))*cos(angle(Y_Bus(n,m))+Delta_in_Rad(m)-Delta_in_Rad(n));
        end
    end
end

% Calculating of J22
J22 = zeros(nbus-PVNum-1,nbus-PVNum-1);
for i=1:nbus-PVNum-1
    n=PQBus(i);
    for jj=1:nbus-PVNum-1
        m=PQBus(jj);
        if n==m
            for k=1:nbus
                J22(i,jj) = J22(i,jj) - V(n)*V(k)*abs(Y_Bus(n,k))*sin(angle(Y_Bus(n,k))+Delta_in_Rad(k)-Delta_in_Rad(n));
            end
            J22(i,jj) = J22(i,jj) - V(n)^2*imag(Y_Bus(n,n));
        else
            J22(i,jj) = J11(n-1,m-1);
        end
    end
end
J=[J11,J12;J21,J22];
end
