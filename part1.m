clc; close all; clear;
syms x1 x2
k = 5;
m = 1;
mi = 0.75;
g = 9.81;
omega = sqrt(k/m);

simulation = [-5 0; 3 0; 0 -5; -6 6; 4 6;];
number = 1; %case simulation 1-5

f1(x1,x2) = (x1 + mi * m * g / k) ^ 2 + x2 ^ 2; 
f2(x1,x2) = (x1 - mi * m * g / k) ^ 2 + x2 ^ 2;

options = odeset('Refine',5);
[t,x]=ode15s(@odefun,[0 5],simulation(number,:),options);

figure(1) %thema 1  
fcontour(f1,[-6 6 0 6])
hold on
fcontour(f2,[-6 6 -6 0]) 
hold on
plot(mi * m * g / k, 0, '+');
hold on
plot(-mi * m * g / k, 0, '+');
title("Phase  Portrait"); xlabel("Position (m)"); ylabel("Velocity (m/s)");

figure(2)
plot(x(:,1),x(:,2))
title("Phase Portrait. Simulation: " + number); xlabel("Position (m)"); ylabel("Velocity (m/s)");

figure(3)
plot(t,x)
title("Simulation: " + number); xlabel("Time (s)"); ylabel("Position (m) / Velocity (m/s)"); legend(["Position","Velocity"]);

function dx=odefun(t,x)
    k = 5;
    m = 1;
    mi = 0.75;
    g = 9.81;
    omega = sqrt(k/m);
    if abs(x(2) * omega) < 10^(-8) % for speed = 0
        if abs(k*x(1)) > abs(m * mi * g) % if the force from the spring is bigger than the static friction
            dx = [0 ;-k*x(1) / (m * omega) - mi * g * sign(x(2) * omega) / omega];
        else % for the area between +-mi*m*g/k
            dx = [0 ;0]; 
        end
            
    else    
        dx=[x(2) * omega ;-k*x(1) / (m * omega) - mi * g * sign(x(2) * omega) / omega];
    end
end