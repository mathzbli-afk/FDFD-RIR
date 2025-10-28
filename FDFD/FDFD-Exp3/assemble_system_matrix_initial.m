function [A,AIP1,AIP2,AIP3] = assemble_system_matrix_initial(kap,h,N,beta)
% -------------------------------------------------------------------------
% Description:
%   Assembles the initial system matrix (A) for the Helmholtz equation
%   using a 27-point finite difference stencil. It separates the matrix
%   into frequency-independent parts (AIP1), parts proportional to k^2
%   (AIP2), and parts proportional to k^4 (AIP3) for the interior points.
%   This allows for efficient recalculation in the frequency sweep.
%
% Inputs:
%   k           - Wavenumber (2*pi*f/c).
%   h           - Grid spacing.
%   N - Vector [Nx, Ny, Nz] of grid cells per dimension.
%   beta        - Vector of 6 boundary admittance coefficients.
%
% Outputs:
%   A           - The complete system matrix for the initial frequency.
%   AIP1, AIP2, AIP3 - Components of the interior points matrix.
% -------------------------------------------------------------------------
NN = N+1;
num_unknowns = NN(1)*NN(2)*NN(3);
A0 = -64/15 + 14*h^2*kap^2/15 - kap^4*h^4/20;
Ass = 7/15 - kap^2*h^2/90;                      % Face neighbors (6)
Asc = 1/10 + kap^2*h^2/90;
Acc = 1/30;
% (2) Predefine neighbor offsets
s_offsets = [1,0,0; -1,0,0; 0,1,0; 0,-1,0; 0,0,1; 0,0,-1];
c_offsets = [1,1,1; 1,1,-1; 1,-1,1; 1,-1,-1; -1,1,1; -1,1,-1; -1,-1,1; -1,-1,-1];
J_offsets = [1,1,0; 1,-1,0; -1,1,0; -1,-1,0; 1,0,1; 1,0,-1; -1,0,1; -1,0,-1; ...
    0,1,1; 0,1,-1; 0,-1,1; 0,-1,-1];
all_offsets = [s_offsets; c_offsets; J_offsets];
all_coeffs = [repmat(Ass, 6, 1); repmat(Acc, 8, 1); repmat(Asc, 12, 1)];


beta_vals_x0 = kap*beta(1);
beta_vals_x1 = kap*beta(2);
beta_vals_y0 = kap*beta(3);
beta_vals_y1 = kap*beta(4);
beta_vals_z0 = kap*beta(5);
beta_vals_z1 = kap*beta(6);


% (4) Pre-allocate sparse matrix triplets (I, J, V)
N_max = 27 * num_unknowns; % Each point has a maximum of 27 connections
% Split into two parts: first, generate the interior points matrix A_IP = A_IP,1 + A_IP,2(k) + A_IP,3(k)
tol = 1;
I1 = zeros(N_max,1);
J1 = zeros(N_max,1);
V_A1 = zeros(N_max,1); V_A2 = zeros(N_max,1); V_A3 = zeros(N_max,1);
for u = 2:NN(1)-1
    for v = 2:NN(2)-1
        for w = 2:NN(3)-1
            % Store the coefficient of the current point
            cnt = cal_cnt(u,v,w,N);
            I1(tol) = cnt; J1(tol) = cnt;
            V_A1(tol) = -64/15; V_A2(tol) = 14/15*kap^2*h^2;V_A3(tol) = -kap^4*h^4/20;
            tol = tol+1;
            % Iterate through the coefficients of neighboring points
            for i=1:6
                offset = s_offsets(i,:);
                u_nei = u + offset(1);
                v_nei = v + offset(2);
                w_nei = w + offset(3);
                col = cal_cnt(u_nei,v_nei,w_nei,N);
                I1(tol) = cnt; J1(tol) = col;
                V_A1(tol) = 7/15; V_A2(tol) = -kap^2*h^2/90; V_A3(tol) = 0;
                tol = tol+1;
            end
            for i=1:8
                offset = c_offsets(i,:);
                u_nei = u + offset(1);
                v_nei = v + offset(2);
                w_nei = w + offset(3);
                col = cal_cnt(u_nei,v_nei,w_nei,N);
                I1(tol) = cnt; J1(tol) = col;
                V_A1(tol) = 1/30; V_A2(tol) = 0; V_A3(tol) = 0;
                tol = tol+1;
            end
            for i=1:12
                offset = J_offsets(i,:);
                u_nei = u + offset(1);
                v_nei = v + offset(2);
                w_nei = w + offset(3);
                col = cal_cnt(u_nei,v_nei,w_nei,N);
                I1(tol) = cnt; J1(tol) = col;
                V_A1(tol) = 1/10; V_A2(tol) = kap^2*h^2/90; V_A3(tol) = 0;
                tol = tol+1;
            end
        end
    end
end
I1 = I1(1:tol-1); J1 = J1(1:tol-1);
V_A1 = V_A1(1:tol-1); V_A2 = V_A2(1:tol-1);V_A3 = V_A3(1:tol-1);
AIP1 = sparse(I1, J1, V_A1, num_unknowns, num_unknowns);
AIP2 = sparse(I1, J1, V_A2, num_unknowns, num_unknowns);
AIP3 = sparse(I1, J1, V_A3, num_unknowns, num_unknowns);
AIP = AIP1+AIP2+AIP3;  % Matrix corresponding to the interior points
% Calculate the matrix corresponding to the boundary points
tol = 1;
I2 = zeros(N_max,1); J2 = zeros(N_max,1);
V_AB = zeros(N_max,1);
for u = 1:NN(1)
    for v = 1:NN(2)
        for w = 1:NN(3)
            if u==1 || u==NN(1) || v==1 || v == NN(2) ||w==1 ||w==NN(3)
                cnt = cal_cnt(u, v, w, N);
                I2(tol) = cnt; J2(tol) = cnt; V_AB(tol) = A0; tol = tol + 1;
                % Iterate through the 26 neighbors
                for s = 1:26
                    offset = all_offsets(s, :);
                    coeff = all_coeffs(s);
                    u_nei = u + offset(1);
                    v_nei = v + offset(2);
                    w_nei = w + offset(3);
                    % Create a task list
                    tasks = [u_nei,v_nei,w_nei,coeff];
                    % Final contribution list for interior points
                    final_contribs = [];
                    while ~isempty(tasks)
                        % Take out the first task
                        current_task = tasks(1,:);
                        tasks(1,:) = [];
                        cu = current_task(1); cv = current_task(2); cw = current_task(3); c_coeff = current_task(4);
                        is_ghost = (cu<1 || cu>NN(1) || cv<1 || cv>NN(2) || cw<1 || cw>NN(3));
                        if ~is_ghost
                            % If it is an interior point, add it to the final contribution list
                            final_contribs = [final_contribs; current_task];
                            continue;
                        end
                        % If it is a ghost point, decompose it and add new tasks to the list
                        if cu<1 % Out of bounds on the left side
                            beta_current = beta_vals_x0;
                            C_term = 2 * 1j * beta_current * h * (1 - (beta_current*h)^2/6 + (beta_current*h)^4/120);
                            tasks = [tasks;[cu+2,cv,cw,c_coeff]];
                            tasks = [tasks;[cu+1,cv,cw,-c_coeff*C_term]];
                        elseif cu>NN(1)
                            beta_current = beta_vals_x1;
                            C_term = 2 * 1j * beta_current * h * (1 - (beta_current*h)^2/6 + (beta_current*h)^4/120);
                            tasks = [tasks;[cu-2,cv,cw,c_coeff]];
                            tasks = [tasks;[cu-1,cv,cw,-c_coeff*C_term]];
                        elseif cv<1
                            beta_current = beta_vals_y0;
                            C_term = 2 * 1j * beta_current * h * (1 - (beta_current*h)^2/6 + (beta_current*h)^4/120);
                            tasks = [tasks;[cu,cv+2,cw,c_coeff]];
                            tasks = [tasks;[cu,cv+1,cw,-c_coeff*C_term]];
                        elseif cv>NN(2)
                            beta_current = beta_vals_y1;
                            C_term = 2 * 1j * beta_current * h * (1 - (beta_current*h)^2/6 + (beta_current*h)^4/120);
                            tasks = [tasks; [cu, cv-2, cw, c_coeff]];
                            tasks = [tasks; [cu, cv-1, cw, -c_coeff * C_term]];
                        elseif cw<1
                            beta_current = beta_vals_z0;
                            C_term = 2 * 1j * beta_current * h * (1 - (beta_current*h)^2/6 + (beta_current*h)^4/120);
                            tasks = [tasks; [cu, cv, cw+2, c_coeff]];
                            tasks = [tasks; [cu, cv, cw+1, -c_coeff * C_term]];
                        elseif cw>NN(3)
                            beta_current = beta_vals_z1;
                            C_term = 2 * 1j * beta_current * h * (1 - (beta_current*h)^2/6 + (beta_current*h)^4/120);
                            tasks = [tasks; [cu, cv, cw-2, c_coeff]];
                            tasks = [tasks; [cu, cv, cw-1, -c_coeff * C_term]];
                        end
                    end
                    % Add all contributions to the sparse matrix
                    for k_contrib = 1:size(final_contribs,1)
                        final_u = final_contribs(k_contrib,1);
                        final_v = final_contribs(k_contrib,2);
                        final_w = final_contribs(k_contrib,3);
                        final_coeff = final_contribs(k_contrib,4);
                        j_col = cal_cnt(final_u,final_v,final_w,N);
                        I2(tol) = cnt; J2(tol) = j_col; V_AB(tol) = final_coeff;
                        tol = tol + 1;
                    end
                end
            end
        end
    end
end
I2 = I2(1:tol-1); J2 = J2(1:tol-1); V_AB = V_AB(1:tol-1);
ABC = sparse(I2,J2,V_AB,num_unknowns,num_unknowns);
A = ABC + AIP;
end