clear;
N = 10^4; % number of symbols
Eb_N0_dB = [0:5:30]; % multiple Eb/N0 values
Nt = 4;
Nr = 4;

ip = (2*(rand(1,N)>0.5)-1) + 1j*(2*(rand(1,N)>0.5)-1);

x_ = reshape(ip, [Nt, N/Nt]);


for Eb_idx = 1:length(Eb_N0_dB)
    disp(Eb_N0_dB(Eb_idx));
    P = sqrt(10^(Eb_N0_dB(Eb_idx)/10))/sqrt(2)/sqrt(Nt);
    x = P * x_; % normalization of energy to P
    cnt_zf = 0;
    cnt_mmse = 0;
    cnt_ml = 0;
    for idx = 1:N/Nt
    % Channel and noise Noise addition
        h = 1/sqrt(2)*[randn(Nr, Nt) + 1j*randn(Nr, Nt)]; % Rayleigh channel
        n = 1/sqrt(2)*[randn(Nr,1) + 1j*randn(Nr,1)];  % white gaussian noise, 0dB variance
        y = h*x(:,idx) + n;


        w_zf = inv(h'* h) * h';
        w_mmse = inv(h'*h + 1/P^2)*h';

        ml_demod = ml_detector(h, y, P, Nt);

        zf_demod_ = qam_demod(w_zf * y);
        zf_demod = reshape(zf_demod_,[Nt, 1]);

        mmse_demod_ = qam_demod(w_mmse * y);
        mmse_demod = reshape(mmse_demod_,[Nt, 1]);


        cnt_ml = cnt_ml + sum(x(:,idx)~=ml_demod,"all");
        cnt_zf = cnt_zf + sum(x_(:,idx)~=zf_demod,"all");
        cnt_mmse = cnt_mmse + sum(x_(:,idx)~=mmse_demod,"all");

    end

    ser_ml(Eb_idx) = cnt_ml/N;
    ser_zf(Eb_idx) = cnt_zf/N;
    ser_mmse(Eb_idx) = cnt_mmse/N;
end

figure
semilogy(Eb_N0_dB, ser_zf, 'v--','Color','#EDB120','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_mmse, '^--','Color','#4DBEEE','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_ml, '-','Color','#000000','LineWidth',2);
legend('ZF', 'MMSE', 'ML');
xlabel('SNR[dB]')
ylabel('SER');
ylim([10^-3.5 10^0]);
title('4 x 4 MIMO, QPSK');


function hat = ml_detector(h, y, P, Nt)
    qam_table = P * [-1-1*1j, 1+1*1j, -1+1*1j, 1-1*1j];
    xx = zeros(Nt, 1);
    x_ = zeros(Nt, 1, length(qam_table).^Nt);

    cnt = 1;
    for idx = 1:length(qam_table)
        xx(1) = qam_table(idx);
        for jdx = 1:length(qam_table)
            xx(2) = qam_table(jdx);
            for kdx = 1:length(qam_table)
                xx(3) = qam_table(kdx);
                for hdx = 1:length(qam_table)
                    xx(4) = qam_table(hdx);
                    x_(:,:,cnt) = xx;
                    result(cnt) = sum(sum(abs(y - h*xx)));
                    cnt = cnt + 1;
                end
            end
        end
    end
    [M I]  = min(result);
    hat = x_(:,:,I);
end


function ipHat = qam_demod(input)
    y_re = real(input);
    y_im = imag(input);
    ipHat(find(y_re < 0 & y_im < 0)) = -1-1*1j;
    ipHat(find(y_re > 0 & y_im > 0)) = 1+1*1j;
    ipHat(find(y_re < 0 & y_im > 0)) = -1+1*1j;
    ipHat(find(y_re > 0 & y_im < 0)) = 1-1*1j;
end