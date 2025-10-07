clc; clear; close all;

% Πραγματικές παράμετροι
m_true = 1.315;
b_true = 0.225;
k_true = 0.725;

% Ρυθμίσεις προσομοίωσης
dt = 0.001;              % μικρό βήμα για ακρίβεια
T = 20;                  % διάρκεια 20 sec
t = 0:dt:T;
N = length(t);

% Είσοδος: ημιτονοειδής
u = 2.5 * sin(t);

% Πραγματικό σύστημα (x, dx, ddx)
x = zeros(1, N); dx = zeros(1, N); ddx = zeros(1, N);
x(1) = 0; dx(1) = 0;
for i = 1:N-1
    ddx(i) = (1/m_true)*(u(i) - b_true*dx(i) - k_true*x(i));
    dx(i+1) = dx(i) + dt * ddx(i);
    x(i+1) = x(i) + dt * dx(i);
end
ddx(end) = ddx(end-1);

% Εκτιμητής (Lyapunov - Parallel Structure)
theta_hat = zeros(3, N);
theta_hat(:,1) = [0.1; 0.5; 0.7];  % αρχικές τιμές κοντά στις σωστές
gamma = 5;                         % ρυθμός προσαρμογής

for i = 1:N-1
    phi = [-dx(i); -x(i); u(i)];
    ddx_hat = phi' * theta_hat(:, i);
    e = ddx(i) - ddx_hat;           % σφάλμα επιτάχυνσης
    theta_hat(:, i+1) = theta_hat(:, i) + gamma * phi * e * dt;
end

% Ανάκτηση φυσικών παραμέτρων
m_hat = 1 ./ theta_hat(3,:);
b_hat = theta_hat(1,:) ./ theta_hat(3,:);
k_hat = theta_hat(2,:) ./ theta_hat(3,:);

% Ανακατασκευή εκτιμημένης κίνησης x̂(t)
x_hat = zeros(1, N); dx_hat = zeros(1, N);
for i = 1:N-1
    ddx_hat_model = (1/m_hat(i)) * (u(i) - b_hat(i)*dx_hat(i) - k_hat(i)*x_hat(i));
    dx_hat(i+1) = dx_hat(i) + dt * ddx_hat_model;
    x_hat(i+1) = x_hat(i) + dt * dx_hat(i+1);
end

% Υπολογισμός σφάλματος εξόδου
e_x = x - x_hat;

% === Γραφήματα ===

% 1. Θέσεις x(t) και x̂(t) και σφάλμα e_x(t)
figure;
subplot(3,1,1);
plot(t, x, 'b', 'LineWidth', 1.5); title('Πραγματική Θέση x(t)');
ylabel('x(t)'); grid on;

subplot(3,1,2);
plot(t, x_hat, 'r', 'LineWidth', 1.5); title('Εκτιμημένη Θέση \hat{x}(t)');
ylabel('\hat{x}(t)'); grid on;

subplot(3,1,3);
plot(t, e_x, 'k', 'LineWidth', 1.5); title('Σφάλμα e_x(t) = x(t) - \hat{x}(t)');
xlabel('Χρόνος (s)'); ylabel('Σφάλμα'); grid on;

% 2. Εκτιμήσεις παραμέτρων
figure;
subplot(3,1,1);
plot(t, m_hat, 'LineWidth', 1.5); yline(m_true, '--r', 'LineWidth', 1.2);
title('Εκτίμηση Μάζας \hat{m}(t)'); ylabel('m(t)'); grid on;

subplot(3,1,2);
plot(t, b_hat, 'LineWidth', 1.5); yline(b_true, '--r', 'LineWidth', 1.2);
title('Εκτίμηση Απόσβεσης \hat{b}(t)'); ylabel('b(t)'); grid on;

subplot(3,1,3);
plot(t, k_hat, 'LineWidth', 1.5); yline(k_true, '--r', 'LineWidth', 1.2);
title('Εκτίμηση Σταθεράς Ελατηρίου \hat{k}(t)'); ylabel('k(t)'); xlabel('Χρόνος (s)'); grid on;
