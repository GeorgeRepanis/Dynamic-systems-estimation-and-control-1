clc;
clear;
close all;

% Πραγματικές παράμετροι συστήματος
a1_real = 1.315;
a2_real = 0.725;
a3_real = 0.225;
b_real = 1.175;

% Παράμετροι εκτιμητή
gamma = diag([0.1, 0.1, 0.1, 0.1]); % Μικρότερο κέρδος για σταθερότητα

% Χρονικό διάστημα
tspan = [0 20];

% Αρχικές συνθήκες: [r(0); r_dot(0); hat_a1(0); hat_a2(0); hat_a3(0); hat_b(0)]
x0 = [0; 0; 0; 0; 0; 0];

% Λύση του συστήματος με ode45
[t, x] = ode45(@(t,x) adaptive_observer_with_disturbance(t, x, a1_real, a2_real, a3_real, b_real, gamma), tspan, x0);

% Εξαγωγή μεταβλητών
r = x(:,1);
r_dot = x(:,2);
hat_a1 = x(:,3);
hat_a2 = x(:,4);
hat_a3 = x(:,5);
hat_b = x(:,6);

% Υπολογισμός επιθυμητής τροχιάς
r_d = (pi/10)*sin(pi/20*t);
r_d_dot = (pi/10)*(pi/20)*cos(pi/20*t);

% Υπολογισμός εκτίμησης r από τις εκτιμήσεις
r_hat = cumtrapz(t, r_dot);

% Υπολογισμός σφάλματος εκτίμησης
e_r = r - r_hat;

% --- Σχεδίαση Γραφημάτων ---

% Γωνίες r(t) και r̂(t)
figure;
plot(t, r, 'b', 'LineWidth', 2);
hold on;
plot(t, r_hat, 'r--', 'LineWidth', 2);
xlabel('Χρόνος [sec]');
ylabel('Γωνία [rad]');
legend('r(t)', 'r̂(t)', 'Location', 'Best');
title('Σύγκριση πραγματικής και εκτιμημένης γωνίας με διαταραχές');
grid on;

% Σφάλμα εκτίμησης
figure;
plot(t, e_r, 'k', 'LineWidth', 2);
xlabel('Χρόνος [sec]');
ylabel('Σφάλμα εκτίμησης e_r(t)');
title('Σφάλμα εκτίμησης γωνίας με διαταραχές');
grid on;

% Εκτιμήσεις παραμέτρων
figure;
plot(t, hat_a1, 'LineWidth', 2);
hold on;
plot(t, hat_a2, 'LineWidth', 2);
plot(t, hat_a3, 'LineWidth', 2);
plot(t, hat_b, 'LineWidth', 2);
xlabel('Χρόνος [sec]');
ylabel('Εκτιμήσεις Παραμέτρων');
legend('hat a_1', 'hat a_2', 'hat a_3', 'hat b', 'Location', 'Best');
title('Εξέλιξη Εκτιμήσεων Παραμέτρων με Διαταραχές');
grid on;

%% --- Συναρτήσεις ---

% Συνάρτηση δυναμικής συστήματος με διαταραχή
function dxdt = adaptive_observer_with_disturbance(t, x, a1, a2, a3, b, gamma)
    r = x(1);
    r_dot = x(2);
    hat_a1 = x(3);
    hat_a2 = x(4);
    hat_a3 = x(5);
    hat_b = x(6);

    % Επιθυμητή τροχιά
    r_d = (pi/10)*sin(pi/20*t);
    r_d_dot = (pi/10)*(pi/20)*cos(pi/20*t);

    % Είσοδος ελέγχου
    u = u_input(t);

    % Εξωτερική διαταραχή
    d = disturbance(t);

    % Πραγματική επιτάχυνση με διαταραχή
    r_ddot = -a1*r_dot - a2*sin(r) + a3*r_dot^2*sin(2*r) + b*u + d;

    % Διάνυσμα φ
    phi = [-r_dot; -sin(r); r_dot^2*sin(2*r); u];

    % Νόμοι προσαρμογής
    theta_hat_dot = gamma * phi * (r - r_d);

    % Διαφορικές εξισώσεις
    dxdt = zeros(6,1);
    dxdt(1) = r_dot;
    dxdt(2) = r_ddot;
    dxdt(3:6) = theta_hat_dot;
end

% Είσοδος ελέγχου
function u = u_input(t)
    u = 0.5 * sin(0.5*t); % Η ίδια είσοδος όπως και πριν
end

% Διαταραχή
function d = disturbance(t)
    d = 0.15 * sin(0.5*t);
end
