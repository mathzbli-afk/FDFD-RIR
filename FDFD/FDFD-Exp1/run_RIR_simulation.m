% run_RIR_simulation.m
% -------------------------------------------------------------------------
% Description:
%   Main script to simulate the Room Impulse Response (RIR) using a
%   sixth-order finite difference method for the Helmholtz equation. This script sets
%   up the simulation parameters, builds the system matrices, iterates
%   through frequencies to solve for the sound pressure, and compares the 
%   frequency-domain RIR results from two different approximated source models.
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
h  = 1/15;            % Grid spacing, uniform in all directions

% -- Room and Transducer Setup
room_dims = [2, 3, 3];   % Room dimensions [Lx, Ly, Lz]
mic_pos   = [0.6, 2.1, 1.3]; % Microphone position [x, y, z]
source_pos= [1.0, 0.5, 1.5]; % Sound source position [x, y, z]

% -- Frequency Range
freq_start = 700;         % Start frequency for simulation
freq_end   = 1000;    % End frequency
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


% Gaussian function and Bump function
for i = 1:2


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
        [b_vector, b1, b2] = generate_source_vector_Gau(k_start,h,N,source_pos);
    else
        [b_vector, b1, b2] = generate_source_vector_Bump(k_start,h,N,source_pos);
    end
    
    % Solve the linear system Ax = b using direct method
    x = A_system \ b_vector;

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
        x = A_system \ b_vector;

        % Reshape and interpolate to get pressure at the microphone
        pressure_field_3D =  reshape(x,N(3)+1,N(2)+1,N(1)+1);
        interpolator = griddedInterpolant({z_nodes,y_nodes,x_nodes},pressure_field_3D,'linear');
        freq_RIR(freq_idx) = interpolator(mic_pos(3), mic_pos(2), mic_pos(1));
    end
    RIR{i}=freq_RIR;
end


TF_Gaus_sliced = RIR{1};
TF_Bump_sliced = RIR{2};



% Create the corresponding frequency vector for the x-axis
frequency_vector = freq_start:freq_end;

% 2. Create the Plot
% Create a new figure window
figure;
hold on; % Allow multiple lines on the same plot

% Plot the Gaussian source data (Solid Blue Line, similar to the reference 'COMSOL' line)
plot(frequency_vector, abs(TF_Gaus_sliced), 'b-', 'LineWidth', 1.5);

% Plot the Bump source data (Dashed Red Line, similar to the reference 'Proposed' line)
plot(frequency_vector, abs(TF_Bump_sliced), 'r--', 'LineWidth', 1.5);

hold off;

% 3. Enhance Aesthetics for Publication
% --- Get the handle to the current axes to customize ---
ax = gca;

% --- Set Title and Labels using LaTeX Interpreter ---
title('Comparison of Approximate Source Models', 'Interpreter', 'latex');
xlabel('Frequency (Hz)', 'Interpreter', 'latex');
ylabel(' Magnitude $|P_h(f)|$', 'Interpreter', 'latex');

% --- Create a Legend ---
legend({'Gaussian Source ($\sigma=h$)', 'Bump Source ($\epsilon=1.5h$)'}, ...
       'Interpreter', 'latex', 'Location', 'northeast');

% --- Customize Grid and Axes ---
grid on; % Turn on the grid
ax.GridLineStyle = '--'; % Set grid lines to be dashed
ax.GridAlpha = 0.5;      % Make grid lines semi-transparent
ax.FontSize = 12;        % Set font size for axis numbers

% --- Set Axis Limits ---
xlim([freq_start, freq_end]);
ylim([0, 0.08]); % Give a little space at the top

% Add a box around the plot area
box on;






