function  plot_rir_results(freq_RIR,fs,p_m,p_s)
% -------------------------------------------------------------------------
% Description:
%   Visualizes the Room Impulse Response (RIR) in both frequency and time
%   domains based on the calculated transfer function.
%
% Inputs:
%   freq_RIR - Complex-valued transfer function at each frequency bin.
%   fs     - The sampling rate (Fs) of the simulation.
%   p_m         - [x, y, z] position of the microphone.
%   p_s      - [x, y, z] position of the sound source.
% -------------------------------------------------------------------------
f_low = 80; % 小于等于80的部分不用处理
freq_RIR(1:f_low)=0;
f_RIR = zeros(fs,1);

for i=2:fs/2+1
    f_RIR(i) = freq_RIR(i-1);
    f_RIR(fs-i+2) = conj(f_RIR(i));
end



t_RIR = real(ifft(f_RIR));
dis_R = sqrt( (p_s(1)-p_m(1))^2 + ...
(p_s(2)-p_m(2))^2 + (p_s(3)-p_m(3))^2 ); % 计算声源和麦克风的距离
coef = ( 1/(4*pi*dis_R) ) *(1/max(t_RIR));
t_RIR = t_RIR*coef;



N = length(t_RIR);
t = (0:N-1) / fs;
f = f_low+1:fs/2;
f_RIR = f_RIR(f_low+2:fs/2+1);



fig = figure('Name', 'RIR Visualization', 'Color', 'white');
fig.Position(3:4) = [700, 600];
subplot(2, 1, 1);
plot(f, abs(f_RIR), 'b-', 'LineWidth', 1.2);
title('Frequency-Domain Room Impulse Response', 'FontSize', 14);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('RIR Magnitude ', 'FontSize', 12);
grid on;
box on;
xlim([f_low, fs/2]);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');
subplot(2, 1, 2);


plot(t, t_RIR, 'r-', 'LineWidth', 1.0);
title('Time-Domain Room Impulse Response', 'FontSize', 14);
xlabel('Time (s)', 'FontSize', 12);
ylabel('RIR Amplitude', 'FontSize', 12);
grid on;
box on;
xlim([0, 1]);
set(gca, 'FontSize', 11, 'FontName', 'Times New Roman');

end