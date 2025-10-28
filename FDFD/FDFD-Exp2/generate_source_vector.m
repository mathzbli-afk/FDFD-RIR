function [b,b1,b2] = generate_source_vector(kap,h,N,p_s,sig)
% -------------------------------------------------------------------------
% Description:
%   Generates the source vector 'b' for the linear system Ax=b.
%   It evaluates the source term and its derivatives (provided as function handles)
%
% Inputs:
%   k              - Wavenumber (2*pi*f/c).
%   h              - Grid spacing.
%   N    - Vector [Nx, Ny, Nz] of grid cells per dimension.
%   p_s, ... - Function handles from create_source_functions.
%
% Outputs:
%   b              - The complete source vector.
%   b1, b2         - Components of b for efficient frequency sweep updates.
% -------------------------------------------------------------------------


% Generate the source vector b
syms x y z  ; % Use symbolic variables to maintain precision
r = (x-p_s(1))^2 + (y-p_s(2))^2 + (z-p_s(3))^2;
F = -4*pi*1/(sig*sqrt(2*pi))*exp(-(r)/(2*sig^2));  % Define the source function, approximated by a Gaussian pulse
Delta_F = laplacian(F,[x,y,z]);
% Sum of the fourth-order pure derivatives
F_four1 = diff(F,x,4)+diff(F,y,4)+diff(F,z,4);
% Sum of the fourth-order mixed derivatives
F_four2 = diff(diff(F,x,2),y,2)+diff(diff(F,x,2),z,2)+diff(diff(F,y,2),z,2);
F_func = matlabFunction(F,'Vars',{x,y,z});
DeltaF_func = matlabFunction(Delta_F,'Vars',{x,y,z});
F_four1_func = matlabFunction(F_four1,'Vars',{x,y,z});
F_four2_func = matlabFunction(F_four2,'Vars',{x,y,z});
% --- Part 2: Efficient computation ---
NN = N + 1; % Number of grid points in each dimension
num_unknowns = NN(1) * NN(2) * NN(3);
b1 = zeros(num_unknowns,1);
b2 = zeros(num_unknowns,1);
id = 1;
for u=1:NN(1)
    for v=1:NN(2)
        for w=1:NN(3)
            X = (u-1)*h;Y=(v-1)*h;Z=(w-1)*h;
            b1(id) = h^2*( F_func(X,Y,Z) + h^2/12*DeltaF_func(X,Y,Z)...
                + h^4/360*F_four1_func(X,Y,Z)+ h^4/90*F_four2_func(X,Y,Z)    );
            b2(id) = -1*h^4*kap^2/20*F_func(X,Y,Z);
            b = b1+b2;
            id = id + 1;
        end
    end
end




end







