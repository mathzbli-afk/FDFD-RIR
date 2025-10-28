function [ABC] = assemble_boundary_matrix(kap, h, N, beta)
% -------------------------------------------------------------------------
% Description:
%   Assembles the boundary condition matrix (ABC) for the Helmholtz
%   equation.
%
% Inputs:
%   k           - Wavenumber (2*pi*f/c).
%   h           - Grid spacing.
%   N - Vector [Nx, Ny, Nz] of grid cells per dimension.
%   beta        - Vector of 6 boundary admittance coefficients.
%
% Output:
%   ABC         - The sparse boundary condition matrix.
% -------------------------------------------------------------------------

NN = N+1; % Number of grid points along each edge
num_unknowns = NN(1)*NN(2)*NN(3);
% --- This part of the code is the same as your original version ---
% (1) Predefine the coefficients for the 27-point stencil
A0 = -64/15 + 14*h^2*kap^2/15 - kap^4*h^4/20;
Ass = 7/15 - kap^2*h^2/90;
Asc = 1/10 + kap^2*h^2/90;
Acc = 1/30;
% (2) Predefine the offsets for the neighbors
s_offsets = [1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
c_offsets = [1,1,1; 1,1,-1; 1,-1,1; 1,-1,-1; -1,1,1; -1,1,-1; -1,-1,1; -1,-1,-1];
J_offsets = [1,1,0; 1,-1,0; -1,1,0; -1,-1,0; 1,0,1; 1,0,-1; -1,0,1; -1,0,-1; ...
    0,1,1; 0,1,-1; 0,-1,1; 0,-1,-1];
all_offsets = [s_offsets; c_offsets; J_offsets];
all_coeffs = [repmat(Ass, 6, 1); repmat(Acc, 8, 1); repmat(Asc, 12, 1)];
% (3) Boundary condition parameters
beta_vals_x0 = kap*beta(1); beta_vals_x1 = kap*beta(2);
beta_vals_y0 = kap*beta(3); beta_vals_y1 = kap*beta(4);
beta_vals_z0 = kap*beta(5); beta_vals_z1 = kap*beta(6);


% =========== Core optimization part: Generate unique boundary point coordinates ===========
% 1. Generate grid point coordinates for each face
[vy, vz] = ndgrid(1:NN(2), 1:NN(3));
face_x1 = [ones(numel(vy), 1), vy(:), vz(:)];         % u = 1
face_xN = [NN(1)*ones(numel(vy), 1), vy(:), vz(:)];    % u = NN(1)
[ux, uz] = ndgrid(1:NN(1), 1:NN(3));
face_y1 = [ux(:), ones(numel(ux), 1), uz(:)];         % v = 1
face_yN = [ux(:), NN(2)*ones(numel(ux), 1), uz(:)];    % v = NN(2)
[ux, uy] = ndgrid(1:NN(1), 1:NN(2));
face_z1 = [ux(:), uy(:), ones(numel(ux), 1)];         % w = 1
face_zN = [ux(:), uy(:), NN(3)*ones(numel(ux), 1)];    % w = NN(3)



% 2. Merge the coordinates of all faces and find the unique coordinate points
all_boundary_coords = [face_x1; face_xN; face_y1; face_yN; face_z1; face_zN];
unique_boundary_coords = unique(all_boundary_coords, 'rows');
num_boundary_points = size(unique_boundary_coords, 1);



% (4) Pre-allocate sparse matrix triplets (I, J, V)
% More accurately estimate the maximum number of non-zero elements
N_max = num_boundary_points * 27;
I2 = zeros(N_max, 1);
J2 = zeros(N_max, 1);
V_AB = zeros(N_max, 1);
tol = 1;
% === New, more efficient loop: Iterate only over unique boundary points ===
for i = 1:num_boundary_points
    u = unique_boundary_coords(i, 1);
    v = unique_boundary_coords(i, 2);
    w = unique_boundary_coords(i, 3);

    % --- The internal processing logic is exactly the same as your original version ---
    cnt = cal_cnt(u, v, w, N);
    I2(tol) = cnt; J2(tol) = cnt; V_AB(tol) = A0; tol = tol + 1;

    % Iterate over the 26 neighbors
    for s = 1:26
        offset = all_offsets(s, :);
        coeff = all_coeffs(s);
        u_nei = u + offset(1);
        v_nei = v + offset(2);
        w_nei = w + offset(3);

        tasks = [u_nei, v_nei, w_nei, coeff];
        final_contribs = [];

        while ~isempty(tasks)
            current_task = tasks(1,:);
            tasks(1,:) = [];
            cu = current_task(1); cv = current_task(2); cw = current_task(3); c_coeff = current_task(4);

            is_ghost = (cu<1 || cu>NN(1) || cv<1 || cv>NN(2) || cw<1 || cw>NN(3));
            if ~is_ghost
                final_contribs = [final_contribs; current_task];
                continue;
            end

            % (The logic for handling ghost points remains unchanged)
            if cu<1 % Exceeded on the left side
                C_term = 2 * 1j * beta_vals_x0 * h * (1 - (beta_vals_x0*h)^2/6 + (beta_vals_x0*h)^4/120);
                tasks = [tasks;[cu+2,cv,cw,c_coeff];[cu+1,cv,cw,-c_coeff*C_term]];
            elseif cu>NN(1)
                C_term = 2 * 1j * beta_vals_x1 * h * (1 - (beta_vals_x1*h)^2/6 + (beta_vals_x1*h)^4/120);
                tasks = [tasks;[cu-2,cv,cw,c_coeff];[cu-1,cv,cw,-c_coeff*C_term]];
            elseif cv<1
                C_term = 2 * 1j * beta_vals_y0 * h * (1 - (beta_vals_y0*h)^2/6 + (beta_vals_y0*h)^4/120);
                tasks = [tasks;[cu,cv+2,cw,c_coeff];[cu,cv+1,cw,-c_coeff*C_term]];
            elseif cv>NN(2)
                C_term = 2 * 1j * beta_vals_y1 * h * (1 - (beta_vals_y1*h)^2/6 + (beta_vals_y1*h)^4/120);
                tasks = [tasks; [cu, cv-2, cw, c_coeff]; [cu, cv-1, cw, -c_coeff * C_term]];
            elseif cw<1
                C_term = 2 * 1j * beta_vals_z0 * h * (1 - (beta_vals_z0*h)^2/6 + (beta_vals_z0*h)^4/120);
                tasks = [tasks; [cu, cv, cw+2, c_coeff]; [cu, cv, cw+1, -c_coeff * C_term]];
            elseif cw>NN(3)
                C_term = 2 * 1j * beta_vals_z1 * h * (1 - (beta_vals_z1*h)^2/6 + (beta_vals_z1*h)^4/120);
                tasks = [tasks; [cu, cv, cw-2, c_coeff]; [cu, cv, cw-1, -c_coeff * C_term]];
            end
        end

        for k_contrib = 1:size(final_contribs,1)
            final_u = final_contribs(k_contrib,1);
            final_v = final_contribs(k_contrib,2);
            final_w = final_contribs(k_contrib,3);
            final_coeff = final_contribs(k_contrib,4);
            j_col = cal_cnt(final_u,final_v,final_w,N);

            if tol > N_max
                I2(end+1000) = 0;
                J2(end+1000) = 0;
                V_AB(end+1000) = 0;
                N_max = length(I2);
            end
            I2(tol) = cnt; J2(tol) = j_col; V_AB(tol) = final_coeff;
            tol = tol + 1;
        end
    end
end
I2 = I2(1:tol-1);
J2 = J2(1:tol-1);
V_AB = V_AB(1:tol-1);
ABC = sparse(I2, J2, V_AB, num_unknowns, num_unknowns);
end