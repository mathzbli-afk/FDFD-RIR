% run_RIR_simulation.m
% -------------------------------------------------------------------------
% Description:
%   Main script to simulate the Room Impulse Response (RIR) using a
%   sixth-order finite difference method for the Helmholtz equation. This script sets
%   up the simulation parameters, builds the system matrices, iterates
%   through frequencies to solve for the sound pressure, and compares the 
%   frequency-domain RIR results under three cases.
%
% Author: Zhibao Li & Wenheng Liu
% Date:   2025.10.12
% -------------------------------------------------------------------------
clear;

% --- 1. Simulation Parameters ---
c_sound = 343.2;      % Speed of sound
rho_0   = 1.2043;     % Density of air

% -- Discretization Parameters
fs = 4000;            % Sampling frequency
h  = 0.1;            % Grid spacing, uniform in all directions

% -- Room and Transducer Setup
room_dims = [2, 3, 3];   % Room dimensions [Lx, Ly, Lz]
mic_pos   = [0.6, 2.1, 1.3]; % Microphone position [x, y, z]
source_pos= [1.0, 0.5, 1.5]; % Sound source position [x, y, z]

% -- Frequency Range
freq_start = 120;         % Start frequency for simulation
freq_end   = 130;    % End frequency
num_freqs  = freq_end - freq_start + 1;

% --- 2. Load Material Impedance Data ---
load("Z_concre_4000.mat"); % Assumed variable name: Z_concre
load("Z_porous_4000.mat"); % Assumed variable name: impedance_Z

% Pre-processing of impedance
low_freq_cutoff_idx = 100;
impedance_Z(1:low_freq_cutoff_idx) = impedance_Z(low_freq_cutoff_idx);
Z_concre(1:low_freq_cutoff_idx)    = Z_concre(low_freq_cutoff_idx);

% --- 3. Grid and System Initialization ---
N = room_dims ./ h;
nodes_per_dim = N + 1;
num_grid_points = prod(nodes_per_dim);

% Create coordinate vectors for each axis
x_nodes = (0:N(1))' * h;
y_nodes = (0:N(2))' * h;
z_nodes = (0:N(3))' * h;

% three cases
for i=1:3

    % Store the frequency domain RIR
    freq_RIR = zeros(num_freqs, 1);

    % --- Calculation for the first frequency (f_start) ---
    f = freq_start;
    k_start = 2 * pi * f / c_sound; % Wavenumber at start frequency

    % Calculate boundary condition coefficients
    beta = calculate_boundary_coeffs(f, impedance_Z, Z_concre);

    % Assemble the full system matrix A and its components for the first frequency
    [A_system, AIP1, AIP2, AIP3] = assemble_system_matrix_initial(k_start, h, N, beta);

    % Generate the source vector b and its components
    if i==1
        [b_vector, b1, b2] = generate_source_vector(k_start,h,N,source_pos,h);
    elseif i==2
        [b_vector, b1, b2] = generate_source_vector(k_start,h,N,source_pos,2*h);
    else
        [b_vector, b1, b2] = generate_source_vector(k_start,h,N,source_pos,h);
    end

    % Solve the linear system Ax = b using GMRES
    if i==1
        x = gmres(A_system, b_vector, 100, 1e-8, 50);
    elseif i==2
        x = gmres(A_system, b_vector, 100, 1e-8, 50);
    else
        x = gmres(A_system, b_vector, 100, 1e-4, 50);
    end


    % Reshape solution vector into a 3D grid and interpolate at mic position
    pressure_field_3D =  reshape(x,N(3)+1,N(2)+1,N(1)+1);
    interpolator = griddedInterpolant({z_nodes,y_nodes,x_nodes},pressure_field_3D,'linear');
    freq_RIR(1) = interpolator(mic_pos(3), mic_pos(2), mic_pos(1));

    % --- Loop for the remaining frequencies ---
    % This loop reuses parts of the matrix A and vector b, which is efficient.
    for freq_idx = 2:num_freqs
        f = freq_start + freq_idx - 1

        k = 2 * pi * f / c_sound; % Current wavenumber

        % Update boundary condition coefficients for the current frequency
        beta = calculate_boundary_coeffs(f, impedance_Z, Z_concre);

        % Assemble only the boundary-dependent part of matrix A
        A_bc = assemble_boundary_matrix(k, h, N, beta);

        % Reconstruct the full system matrix using pre-calculated components
        % This avoids re-calculating the entire interior matrix every time.
        k_ratio = (k / k_start)^2;
        A_system = A_bc + AIP1 + k_ratio * AIP2 + (k_ratio^2) * AIP3;

        % Reconstruct the source vector
        b_vector = b1 + k_ratio * b2; % Correction based on original code structure

        % Solve the linear system
        if i==1
            x = gmres(A_system, b_vector, 100, 1e-8, 50);
        elseif i==2
            x = gmres(A_system, b_vector, 100, 1e-8, 50);
        else
            x = gmres(A_system, b_vector, 100, 1e-4, 50);
        end

        % Reshape and interpolate to get pressure at the microphone
        pressure_field_3D =  reshape(x,N(3)+1,N(2)+1,N(1)+1);
        interpolator = griddedInterpolant({z_nodes,y_nodes,x_nodes},pressure_field_3D,'linear');
        freq_RIR(freq_idx) = interpolator(mic_pos(3), mic_pos(2), mic_pos(1));
    end
    RIR{i} = freq_RIR;
end


%--- 5. Post-Processing and Visualization ---


freq = freq_start:freq_end;
mag_matlab = abs(RIR{1})./max(abs(RIR{1}));
mag_matlab_2h = abs(RIR{2})./max(abs(RIR{2}));
mag_matlab_h_tol_4 = abs(RIR{3})./max(abs(RIR{3}));


% 2. Create figure and plot
% Create a new figure
figure('Position', [100, 100, 800, 450]);
hold on; % Allow multiple plots on the same axes
% Plot the first algorithm curve (Case 1: red, dashed line)
plot(freq, mag_matlab, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Case 1 ($\sigma=h$, $\epsilon_{tol}=10^{-8}$)');
% Plot the second algorithm curve (Case 2: dark green, dash-dot line)
dark_green = [0, 0.6, 0.2];
plot(freq, mag_matlab_2h, '-.', 'Color', dark_green, 'LineWidth', 1.5, 'DisplayName', 'Case 2 ($\sigma=2h, \epsilon_{tol}=10^{-8}$)');
% --- Added: Plot the third curve (Case 3) ---
plot(freq, mag_matlab_h_tol_4, 'm:', 'LineWidth', 2, 'DisplayName', 'Case 3 ($\sigma=h$, $\epsilon_{tol}=10^{-4}$)');
% --- End of addition ---
hold off;


% 3. Add formatting for publication
% Set title and axis labels, using LaTeX interpreter
title('Magnitude Comparison', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('Frequency (Hz)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Normalized Magnitude $|P_h(f)|$', 'Interpreter', 'latex', 'FontSize', 14);
% Create and set the legend
lgd = legend('show');
set(lgd, 'Location', 'NorthEast', 'Interpreter', 'latex', 'FontSize', 12);
% Set axis tick font and font size
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
% Add grid and border
grid on;
box on;
set(gca, 'GridAlpha', 0.4, 'GridLineStyle', '--');
% Set the X-axis limits
xlim([min(freq), max(freq)]);


% --- Add inset (zoomed-in plot) ---
% Create a small axes, positioned in the bottom-right of the main plot
inset_pos = [0.6, 0.4, 0.25, 0.3 ];
axes('Position', inset_pos); % [left, bottom, width, height]
hold on;
plot(freq, mag_matlab, 'r--', 'LineWidth', 1.5);
plot(freq, mag_matlab_2h, '-.', 'Color', dark_green, 'LineWidth', 1.5);
plot(freq, mag_matlab_h_tol_4, 'm:', 'LineWidth', 2);
hold off;
box on;
% Set the axis limits for the zoomed-in region, e.g., zoom in on the 110-130Hz peak
xlim([110, 130]);
ylim([0.4, 0.6]); % Adjust according to your data
set(gca, 'FontSize', 10);
annotation('rectangle', [0.4, 0.43, 0.05, 0.06], 'Color', 'k', 'LineStyle', '--');