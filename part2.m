clc; close all; clear;

options = odeset('Refine',5);
[t,x]=ode15s(@odefun,[0 15],[3 0] ,options);
[t2,x2]=ode15s(@odefun2,[0 15],[3 0] ,options);
[t3,x3]=ode15s(@odefun3,[0 15],[3 0] ,options);
[t4,x4]=ode15s(@odefun4,[0 15],[3 0] ,options);

figure(1)
plot(x(:,1),x(:,2))
title("Phase  Portrait. Sliding Mode Control"); xlabel("Position (m)"); ylabel("Velocity (m/s)");
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
axis equal
figure(2)
plot(t,x)
title("Sliding Mode Control"); xlabel("Time (s)"); ylabel("Position (m) / Velocity (m/s)"); legend(["Position","Velocity"]);

figure(3)
plot(x2(:,1),x2(:,2))
title("Phase  Portrait. Sliding Mode Control without f(x)"); xlabel("Position (m)"); ylabel("Velocity (m/s)");

figure(4)
plot(t2,x2)
title("Sliding Mode Control without f(x)"); xlabel("Time (s)"); ylabel("Position (m) / Velocity (m/s)"); legend(["Position","Velocity"]);

figure(5)
plot(x3(:,1),x3(:,2))
title("Phase  Portrait. Lyapunov Redesign"); xlabel("Position (m)"); ylabel("Velocity (m/s)");

figure(6)
plot(t3,x3)
title("Lyapunov Redesign"); xlabel("Time (s)"); ylabel("Position (m) / Velocity (m/s)"); legend(["Position","Velocity"]);

figure(7)
plot(x4(:,1),x4(:,2))
title("Phase  Portrait. Lyapunov Redesign without f(x)"); xlabel("Position (m)"); ylabel("Velocity (m/s)");

figure(8)
plot(t4,x4)
title("Lyapunov Redesign without f(x)"); xlabel("Time (s)"); ylabel("Position (m) / Velocity (m/s)"); legend(["Position","Velocity"]);


function dx=odefun(t,x)
    m = 1;
    mi = 0.75;
    g = 9.81;
    lambda = 1;
    c = 1;
    xd = 2 + sin(2.5 * t) + 2 * cos(1.25 * t); 
    dxd = 2.5 * cos(2.5 * t) - 2.5 * sin(1.25 * t);
    ddxd = - 6.25 * sin(2.5 * t) -3.125 * cos(1.25 * t);
    e = x(1) - xd;
    de = x(2)  - dxd;
    p = 0.9375 * g + 3.5 * abs(x(1)) + 4.5 * abs(x(1))^3 + 0.75 * abs(ddxd - lambda * de) + c;
    f = (5 + 0.6 * exp(-0.045 * t)) * x(1) * (1 + 0.9 * x(1)^2);
    u = 1.0625 * g * sign(x(2)  ) + 7.5 * x(1) + 6.5 * x(1)^3 + 1.25*(ddxd - lambda * de) - p * sat(de + lambda * e);
    if abs(x(2) ) < 10^(-8)
        if abs(u - f) > abs(m * mi * g)
            dx = [0;  1/m * (u - m * mi * g * sign(x(2)) - f)];
        else
            dx = [0;0];
        end
    else
        dx = [x(2) ; 1/m * (u - m * mi * g * sign(x(2)) - f)];
    end
end

function dx=odefun2(t,x)
    m = 1;
    mi = 0.75;
    g = 9.81;
    lambda = 1;
    c = 1;
    xd = 2 + sin(2.5 * t) + 2 * cos(1.25 * t); 
    dxd = 2.5 * cos(2.5 * t) - 2.5 * sin(1.25 * t);
    ddxd = - 6.25 * sin(2.5 * t) -3.125 * cos(1.25 * t);
    e = x(1) - xd;
    de = x(2)  - dxd;
    p = 0.9375 * g + 22 * abs(x(1))^3 + 0.75 * abs(ddxd - lambda * de) + c;
    f = (5 + 0.6 * exp(-0.045 * t)) * x(1) * (1 + 0.9 * x(1)^2);
    u = 1.0625 * g * sign(x(2) ) + 1.25*(ddxd - lambda * de) - p * sat(de + lambda * e);
    if abs(x(2) ) < 10^(-8)
        if abs(u - f) > abs(m * mi * g)
            dx = [0;  1/m * (u - m * mi * g * sign(x(2)) - f)];
        else
            dx = [0;0];
        end
    else
        dx = [x(2) ; 1/m * (u - m * mi * g * sign(x(2)) - f)];
    end
end

function dx=odefun3(t,x)
    m = 1;
    mi = 0.75;
    g = 9.81;
    xd = 2 + sin(2.5 * t) + 2 * cos(1.25 * t); 
    dxd = 2.5 * cos(2.5 * t) - 2.5 * sin(1.25 * t);
    f = (5 + 0.6 * exp(-0.045 * t)) * x(1) * (1 + 0.9 * x(1)^2);
    e = [x(1) - xd;x(2) - dxd];
    F = [1;1];
    P = [3/2  1/2 ; 1/2  1];
    b = [0;1];
    v = -transpose(F) * e;
    p = 1/0.375 *( 0.625 * abs(v) + 2 * (1.6484375 * g + 3.5 * abs(x(1)) + 4.5 * abs(x(1)) ^ 3) );
    dv = -p  * sat(transpose(e) * P * b);
    u = 0.3515625 * g * sign(x(2)) + 7.5 * x(1) + 6.5 * x(1)^3 + (v + dv) * 0.75;
    if abs(x(2) ) < 10^(-8)
        if abs(u - f) > abs(m * mi * g)
            dx = [0;  1/m * (u - m * mi * g * sign(x(2)) - f)];
        else
            dx = [0;0];
        end
    else
        dx = [x(2); 1/m * (u - m * mi * g * sign(x(2)) - f)];
    end
end

function dx=odefun4(t,x)
    m = 1;
    mi = 0.75;
    g = 9.81;
    xd = 2 + sin(2.5 * t) + 2 * cos(1.25 * t); 
    dxd = 2.5 * cos(2.5 * t) - 2.5 * sin(1.25 * t);
    f = (5 + 0.6 * exp(-0.045 * t)) * x(1) * (1 + 0.9 * x(1)^2);
    e = [x(1) - xd;x(2) - dxd];
    F = [1;1];
    P = [3/2  1/2 ; 1/2  1];
    b = [0;1];
    v = -transpose(F) * e;
    p = 1/0.375 *(0.625 * abs(v) + 2 * (1.6484375 * g + 22 * abs(x(1)) ^ 3) );
    dv = -p * sat(transpose(e) * P * b);
    u = 0.3515625 * g * sign(x(2)) + (v + dv) * 0.75;
    if abs(x(2) ) < 10^(-8)
        if abs(u - f) > abs(m * mi * g)     
            dx = [0;  1/m * (u - m * mi * g * sign(x(2)) - f)];
        else
            dx = [0;0];
        end
    else
        dx = [x(2); 1/m * (u - m * mi * g * sign(x(2)) - f)];
    end
end

function y = sat(s)
    if s > 1
        y = 1;
    elseif s < -1 
        y = -1;
    else 
        y = s;
    end
end