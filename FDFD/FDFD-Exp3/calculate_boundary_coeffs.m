function beta = calculate_boundary_coeffs(f,impedance_Z,Z_concre)
    beta = zeros(6,1);

    Z_p = impedance_Z(f);
    Z_c = Z_concre(f);
    beta(1) = Z_p;
    beta(2:6) = Z_c;
end
