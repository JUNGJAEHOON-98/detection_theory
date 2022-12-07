clear;
N = 10^6; % number of symbols
Eb_N0_dB = [0:5:30]; % multiple Eb/N0 values
Nt = 4;
Nr = 4;

ip = (2*(rand(1,N)>0.5)-1) + 1j*(2*(rand(1,N)>0.5)-1);
s = (1/sqrt(2))*ip; % normalization of energy to 1
x = reshape(s,[Nt, N/Nt]);

h = 1/sqrt(2) * [randn(Nr, Nt) + 1j*randn(Nr, Nt)]; % Rayleigh channel
n = 1/sqrt(2) * [randn(Nr, N/Nt) + 1j*randn(Nr, N/Nt)]; % white gaussian noise, 0dB variance

for idx = 1:length(Eb_N0_dB)

    % Channel and noise Noise addition
    y = h*x + 10^(-Eb_N0_dB(idx)/20)*n ;
    w_zf = inv(h'*h)*h';
    w_mmse = inv(h'*h + sqrt(10^(-Eb_N0_dB(idx)/10))/Nt * eye(Nt))*h';

    zf_demod = qam_demod(w_zf * y);
    mmse_demod = qam_demod(w_mmse * y);

    ser_zf(idx) = sum(ip~=zf_demod,"all")/N;
    ser_mmse(idx) = sum(ip~=mmse_demod,"all")/N;
end

figure
semilogy(Eb_N0_dB, ser_zf, 'bp-');
hold on
semilogy(Eb_N0_dB, ser_mmse, 'mo-');
legend('ZF', 'MMSE');
xlabel('SNR[dB]')
ylabel('SER');
title('4 x 4 MIMO, QPSK');


function ipHat = qam_demod(input)
    y_re = real(input);
    y_im = imag(input);
    ipHat(find(y_re < 0 & y_im < 0)) = -1 + -1*1j;
    ipHat(find(y_re >= 0 & y_im > 0)) = 1 + 1*1j;
    ipHat(find(y_re < 0 & y_im >= 0)) = -1 + 1*1j;
    ipHat(find(y_re >= 0 & y_im < 0)) = 1 - 1*1j;
end