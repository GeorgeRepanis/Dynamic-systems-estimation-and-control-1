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

% Είσοδος του συστήματος
u = 2.5 * sin(t);

% Διαφορετικές τιμές πλάτους θορύβου
eta0_values = 0:0.05:0.5;   
num_cases = length(eta0_values);

% Πίνακες για αποθήκευση σφαλμάτων
error_m = zeros(1, num_cases);
error_b = zeros(1, num_cases);
error_k = zeros(1, num_cases);

% === Βρόχος για κάθε πλάτος θορύβου ===
for j = 1:num_cases
    eta0 = eta0_values(j);
    f0 = 20;  % Συχνότητα θορύβου
    eta = eta0 * sin(2*pi*f0*t);  % Θόρυβος

    % --- Πραγματικό σύστημα ---
    x = zeros(1, N); dx = zeros(1, N); ddx = zeros(1, N);
    x(1) = 0; dx(1) = 0;
    for i = 1:N-1
        ddx(i) = (1/m_true)*(u(i) - b_true*dx(i) - k_true*x(i));
        dx(i+1) = dx(i) + dt * ddx(i);
        x(i+1) = x(i) + dt * dx(i);
    end
    ddx(end) = ddx(end-1);

    % --- Θορυβώδη μέτρηση ---
    x_noise = x + eta;
    dx_noise = gradient(x_noise, dt);
    ddx_noise = gradient(dx_noise, dt);

    % --- Εκτιμητής Μεικτής Δομής ---
    theta_hat = zeros(3, N);
    theta_hat(:,1) = [0.1; 0.5; 0.7];  % Αρχικές εκτιμήσεις
    gamma = 5;
    beta = 2;
    
    x_hat = zeros(1, N); dx_hat = zeros(1, N);

    for i = 1:N-1
        phi = [-dx_noise(i); -x_noise(i); u(i)];
        ddx_hat = phi' * theta_hat(:, i);
        e_ddx = ddx_noise(i) - ddx_hat;
        e_x = x(i) - x_hat(i);

        theta_hat(:, i+1) = theta_hat(:, i) + gamma * (phi * e_ddx + beta * phi * e_x) * dt;
        
        % Ανακατασκευή κίνησης
        m_hat_i = 1/theta_hat(3,i);
        b_hat_i = theta_hat(1,i)/theta_hat(3,i);
        k_hat_i = theta_hat(2,i)/theta_hat(3,i);

        ddx_hat_model = (1/m_hat_i)*(u(i) - b_hat_i*dx_hat(i) - k_hat_i*x_hat(i));
        dx_hat(i+1) = dx_hat(i) + dt*ddx_hat_model;
        x_hat(i+1) = x_hat(i) + dt*dx_hat(i+1);
    end

    % --- Ανάκτηση φυσικών παραμέτρων ---
    m_hat = 1 ./ theta_hat(3,:);
    b_hat = theta_hat(1,:) ./ theta_hat(3,:);
    k_hat = theta_hat(2,:) ./ theta_hat(3,:);

    % --- Υπολογισμός σφαλμάτων στο τέλος της προσομοίωσης ---
    error_m(j) = abs(m_true - m_hat(end));
    error_b(j) = abs(b_true - b_hat(end));
    error_k(j) = abs(k_true - k_hat(end));
end

% === Διάγραμμα Σφαλμάτων ===

figure;
plot(eta0_values, error_m, '-o', 'LineWidth', 1.5); hold on;
plot(eta0_values, error_b, '-s', 'LineWidth', 1.5);
plot(eta0_values, error_k, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Πλάτος Θορύβου \eta_0');
ylabel('Σφάλμα Εκτίμησης');
title('Εξέλιξη Σφαλμάτων Εκτίμησης Παραμέτρων σε Συνάρτηση με το Πλάτος Θορύβου (Μεικτή Δομή)');
legend('Σφάλμα Μάζας m', 'Σφάλμα Απόσβεσης b', 'Σφάλμα Σταθεράς k');
