clc;
clear;

% Input Variables
L = 0.3;
rho = 2500;
c = 840;
k = 40;
alpha = 1 / (rho * c);
h = 25; % Convective heat transfer coefficient

% Grid Generation
N = 100;
deltax = L / N;
x = (deltax / 2):deltax:L;
x1 = 0:deltax:L;

% Critical Time Determination
delta_t_critical = ((1 / alpha) * (deltax ^ 2)) / (2 * k);

% Time Step Input
T = input('Please enter the time (s): ', 's');
t = T*200;
fprintf('The critical time step value is %f\n', delta_t_critical);
deltat = 0.1;
M = (t / deltat) + 1;
time = zeros(1, M);
for a = 1:M
    time(a) = (a - 1) * deltat;
end

% Initial Conditions
tini1 = 6;
tini2 = 6;
temp = zeros(N, M);
for i = 1:N
    j = 1;
    if x(i) <= (L / 2)
        temp(i, j) = tini1;
    else
        temp(i, j) = tini1 + (((x(i) - L / 2) * (tini2 - tini1)) / (L - (L / 2)));
    end
end

% Boundary Conditions
qin = 333333; % Heat flux at x = L

% Solving
for j = 2:M
    for i = 1:N
        if x(i) == (deltax / 2)
            temp(i, j) = alpha * (deltat / deltax) * ((k * temp(i + 1, j - 1) / (x(i + 1) - x(i))) + (k * temp(i, j - 1) / (x(i + 1) - x(i))) + (((inv(alpha) * deltax / deltat) - (k / (x(i + 1) - x(i))) - (k / (x(i + 1) - x(i)))) * temp(i, j - 1))) + 2*(h * (temp(i, j - 1) - 14) * deltat / (rho * c));
        elseif x(i) == (L - (deltax / 2))
            temp(i, j) = alpha * (deltat / deltax) * (((k * temp(i, j - 1) + (qin * (deltax / 2) / k)) / (L - x(i))) + (k * temp(i - 1, j - 1) / (x(i) - x(i - 1))) + (((inv(alpha) * deltax / deltat) - (k / (x(i) - x(i - 1))) - (k / (L - x(i)))) * temp(i, j - 1))) + 2*(h * (temp(i, j - 1) - 14) * deltat / (rho * c));
        else
            temp(i, j) = alpha * (deltat / deltax) * ((k * temp(i + 1, j - 1) / (x(i + 1) - x(i))) + (k * temp(i - 1, j - 1) / (x(i) - x(i - 1))) + (((inv(alpha) * deltax / deltat) - (k / (x(i + 1) - x(i))) - (k / (x(i) - x(i - 1)))) * temp(i, j - 1))) + 2*(h * (temp(i, j - 1) - 14) * deltat / (rho * c));
        end
    end
end

% Post Processing
plot(time, temp);
title('Temperature vs Time');
grid on;
xlabel('Time in seconds');
ylabel('Temperature in Celcius');

figure;
plot(x, temp(:, end))
grid on;
title('Temperature vs Length');
xlabel('Length in meters');
ylabel('Temperature in Celcius');
