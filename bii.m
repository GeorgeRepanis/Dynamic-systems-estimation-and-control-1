clc; clear; close all;

% Πραγματικές παράμετροι
m_true = 1.315;
b_true = 0.225;
k_true = 0.725;

% Ρυθμίσεις προσομοίωσης
dt = 0.001;           
T = 20;               
t = 0:dt:T;
N = length(t);

% Είσοδος
u = 2.5 * sin(t);

% === Πραγματικό σύστημα χωρίς θόρυβο ===
x = zeros(1, N); dx = zeros(1, N); ddx = zeros(1, N);
x(1) = 0; dx(1) = 0;
for i = 1:N-1
    ddx(i) = (1/m_true)*(u(i) - b_true*dx(i) - k_true*x(i));
    dx(i+1) = dx(i) + dt * ddx(i);
    x(i+1) = x(i) + dt * dx(i);
end
ddx(end) = ddx(end-1);

% Προσθήκη θορύβου στην έξοδο
eta0 = 0.25;   
f0 = 20;      
eta = eta0 * sin(2*pi*f0*t); 
x_noise = x + eta;  % Θορυβώδες x

% === Εκτιμητής Παράλληλης Δομής χωρίς θόρυβο ===
theta_hat_clean = zeros(3, N);
theta_hat_clean(:,1) = [0.1; 0.5; 0.7];
gamma = 5;

x_hat_clean = zeros(1, N); dx_hat_clean = zeros(1, N);
for i = 1:N-1
    phi = [-dx(i); -x(i); u(i)];
    ddx_hat = phi' * theta_hat_clean(:, i);
    e = ddx(i) - ddx_hat;
    theta_hat_clean(:, i+1) = theta_hat_clean(:, i) + gamma * phi * e * dt;
    
    % Ανακατασκευή κίνησης
    m_hat_i = 1/theta_hat_clean(3,i);
    b_hat_i = theta_hat_clean(1,i)/theta_hat_clean(3,i);
    k_hat_i = theta_hat_clean(2,i)/theta_hat_clean(3,i);
    
    ddx_hat_model = (1/m_hat_i)*(u(i) - b_hat_i*dx_hat_clean(i) - k_hat_i*x_hat_clean(i));
    dx_hat_clean(i+1) = dx_hat_clean(i) + dt*ddx_hat_model;
    x_hat_clean(i+1) = x_hat_clean(i) + dt*dx_hat_clean(i+1);
end

% === Εκτιμητής Παράλληλης Δομής με θόρυβο ===
theta_hat_noise = zeros(3, N);
theta_hat_noise(:,1) = [0.1; 0.5; 0.7];

x_hat_noise = zeros(1, N); dx_hat_noise = zeros(1, N);
dx_noise = gradient(x_noise, dt);
ddx_noise = gradient(dx_noise, dt);

for i = 1:N-1
    phi_noise = [-dx_noise(i); -x_noise(i); u(i)];
    ddx_hat_noise = phi_noise' * theta_hat_noise(:, i);
    e_noise = ddx_noise(i) - ddx_hat_noise;
    theta_hat_noise(:, i+1) = theta_hat_noise(:, i) + gamma * phi_noise * e_noise * dt;
    
    % Ανακατασκευή κίνησης
    m_hat_i = 1/theta_hat_noise(3,i);
    b_hat_i = theta_hat_noise(1,i)/theta_hat_noise(3,i);
    k_hat_i = theta_hat_noise(2,i)/theta_hat_noise(3,i);
    
    ddx_hat_model_noise = (1/m_hat_i)*(u(i) - b_hat_i*dx_hat_noise(i) - k_hat_i*x_hat_noise(i));
    dx_hat_noise(i+1) = dx_hat_noise(i) + dt*ddx_hat_model_noise;
    x_hat_noise(i+1) = x_hat_noise(i) + dt*dx_hat_noise(i+1);
end

% Ανάκτηση φυσικών παραμέτρων
m_hat_clean = 1 ./ theta_hat_clean(3,:);
b_hat_clean = theta_hat_clean(1,:) ./ theta_hat_clean(3,:);
k_hat_clean = theta_hat_clean(2,:) ./ theta_hat_clean(3,:);

m_hat_noise = 1 ./ theta_hat_noise(3,:);
b_hat_noise = theta_hat_noise(1,:) ./ theta_hat_noise(3,:);
k_hat_noise = theta_hat_noise(2,:) ./ theta_hat_noise(3,:);

% Σφάλματα εξόδου
e_x_clean = x - x_hat_clean;
e_x_noise = x - x_hat_noise;

% === Γραφήματα ===

% Θέσεις και σφάλματα
figure;
subplot(3,1,1);
plot(t, x, 'b', t, x_hat_clean, 'r--', 'LineWidth', 1.5);
title('Χωρίς Θόρυβο: x(t) και \hat{x}(t)');
ylabel('Θέση'); grid on; legend('x(t)','\hat{x}(t)');

subplot(3,1,2);
plot(t, x, 'b', t, x_hat_noise, 'r--', 'LineWidth', 1.5);
title('Με Θόρυβο: x(t) και \hat{x}(t)');
ylabel('Θέση'); grid on; legend('x(t)','\hat{x}(t)');

subplot(3,1,3);
plot(t, e_x_clean, 'b', t, e_x_noise, 'r--', 'LineWidth', 1.5);
title('Σφάλμα Εξόδου: e_x(t)');
xlabel('Χρόνος (s)'); ylabel('Σφάλμα'); grid on; legend('Χωρίς Θόρυβο','Με Θόρυβο');

% Εκτίμηση Παραμέτρων
figure;
subplot(3,1,1);
plot(t, m_hat_clean, 'b', t, m_hat_noise, 'r--', 'LineWidth', 1.5); hold on;
yline(m_true, '--k', 'LineWidth', 1.2);
title('Εκτίμηση Μάζας \hat{m}(t)');
ylabel('m(t)'); grid on; legend('Χωρίς Θόρυβο','Με Θόρυβο','Πραγματική');

subplot(3,1,2);
plot(t, b_hat_clean, 'b', t, b_hat_noise, 'r--', 'LineWidth', 1.5); hold on;
yline(b_true, '--k', 'LineWidth', 1.2);
title('Εκτίμηση Απόσβεσης \hat{b}(t)');
ylabel('b(t)'); grid on; legend('Χωρίς Θόρυβο','Με Θόρυβο','Πραγματική');

subplot(3,1,3);
plot(t, k_hat_clean, 'b', t, k_hat_noise, 'r--', 'LineWidth', 1.5); hold on;
yline(k_true, '--k', 'LineWidth', 1.2);
title('Εκτίμηση Σταθεράς Ελατηρίου \hat{k}(t)');
xlabel('Χρόνος (s)'); ylabel('k(t)'); grid on; legend('Χωρίς Θόρυβο','Με Θόρυβο','Πραγματική');
