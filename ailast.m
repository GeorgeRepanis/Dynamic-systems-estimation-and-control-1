clc; clear; close all;

% Πραγματικές τιμές παραμέτρων
m_true = 1.315;
b_true = 0.225;
k_true = 0.725;

% Ρυθμίσεις προσομοίωσης
dt = 0.001;
T = 20;
t = 0:dt:T;
N = length(t);

% Είσοδος: u(t) = 2.5
u = 2.5 * ones(1, N);

% Πραγματική κίνηση x(t)
x = zeros(1, N); dx = zeros(1, N); ddx = zeros(1, N);
x(1) = 0; dx(1) = 0;
for i = 1:N-1
    ddx(i) = (1/m_true)*(u(i) - b_true*dx(i) - k_true*x(i));
    dx(i+1) = dx(i) + dt * ddx(i);
    x(i+1) = x(i) + dt * dx(i);
end
ddx(end) = ddx(end-1);

% Εκτιμητής - Gradient
theta_hat = zeros(3, N);
theta_hat(:,1) = [0.1; 0.5; 0.7];  % Αρχικές εκτιμήσεις
gamma = 0.5;                       % Ρυθμός μάθησης

for i = 1:N-1
    phi = [-dx(i); -x(i); u(i)];
    ddx_hat = phi' * theta_hat(:, i);
    e = ddx(i) - ddx_hat;
    
    theta_hat(:, i+1) = theta_hat(:, i) + gamma * phi * e * dt;
    
    % Προστασία θ3 για αποφυγή αριθμητικής κατάρρευσης
    if abs(theta_hat(3,i+1)) < 0.05
        theta_hat(3,i+1) = 0.05;
    end
end

% Υπολογισμός παραμέτρων m, b, k
m_hat = 1 ./ theta_hat(3,:);
b_hat = theta_hat(1,:) ./ theta_hat(3,:);
k_hat = theta_hat(2,:) ./ theta_hat(3,:);

% Περιορισμός (saturation) στις παραμέτρους
m_hat = min(max(m_hat, 1), 2);       % [1, 2]
b_hat = min(max(b_hat, 0), 0.5);     % [0, 0.5]
k_hat = min(max(k_hat, 0.5), 1);     % [0.5, 1]

% Ανακατασκευή εκτιμώμενης κίνησης x̂(t)
x_hat = zeros(1, N);
dx_hat = zeros(1, N);
ddx_hat = zeros(1, N);
for i = 1:N-1
    ddx_hat(i) = (1/m_hat(i)) * (u(i) - b_hat(i)*dx_hat(i) - k_hat(i)*x_hat(i));
    dx_hat(i+1) = dx_hat(i) + dt * ddx_hat(i);
    x_hat(i+1) = x_hat(i) + dt * dx_hat(i);
end
ddx_hat(end) = ddx_hat(end-1);
e_x = x - x_hat;

% === Γραφήματα ===
figure;
subplot(3,1,1); plot(t, x, 'b', 'LineWidth', 1.5); title('x(t)');
subplot(3,1,2); plot(t, x_hat, 'r', 'LineWidth', 1.5); title('hat{x}(t)');
subplot(3,1,3); plot(t, e_x, 'k', 'LineWidth', 1.5); title('e_x(t) = x(t) - hat{x}(t)');

figure;
subplot(3,1,1); plot(t, m_hat, 'LineWidth', 1.5); yline(m_true, '--r');
title('Εκτίμηση μάζας hat/{m}(t)'); ylabel('m');

subplot(3,1,2); plot(t, b_hat, 'LineWidth', 1.5); yline(b_true, '--r');
title('Εκτίμηση απόσβεσης hat/{b}(t)'); ylabel('b');

subplot(3,1,3); plot(t, k_hat, 'LineWidth', 1.5); yline(k_true, '--r');
title('Εκτίμηση σταθεράς ελατηρίου hat/{k}(t)'); ylabel('k'); xlabel('Χρόνος (s)');

