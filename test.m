clear;
P = 1;
N = 10^5; % number of symbols
Eb_N0_dB = [0:6:30]; % multiple Eb/N0 values
Nt = 4;
Nr = 4;

ip = (2*(rand(1,N)>0.5)-1) + 1j*(2*(rand(1,N)>0.5)-1);

x_ = reshape(ip, [Nt, N/Nt]);
x = P * (1/sqrt(2))*x_; % normalization of energy to P

h = 1/sqrt(2)*[randn(Nr, Nt, N/Nt) + 1j*randn(Nr, Nt, N/Nt)]; % Rayleigh channel
n = 1/sqrt(2)*[randn(Nr, N/Nt) + 1j*randn(Nr, N/Nt)]; % white gaussian noise, 0dB variance

for Eb_idx = 1:length(Eb_N0_dB)
    disp(Eb_N0_dB(Eb_idx));
    cnt_zf = 0;
    cnt_mmse = 0;
    cnt_ml = 0;
    for idx = 1:N/Nt
    % Channel and noise Noise addition

        y = h(:,:,idx)*x(:,idx) + 10^(-Eb_N0_dB(Eb_idx)/20)*n(:,idx);

        w_zf = inv(h(:,:,idx)'* h(:,:,idx)) * h(:,:,idx)';
        w_mmse = inv(h(:,:,idx)'*h(:,:,idx) + 10^(-Eb_N0_dB(Eb_idx)/10)/Nt * eye(Nt))*h(:,:,idx)';

        ml_demod = ml_detector(h(:,:,idx), y, Nt);

        zf_demod_ = qam_demod(w_zf * y);
        zf_demod = reshape(zf_demod_,[Nt, 1]);

        mmse_demod_ = qam_demod(w_mmse * y);
        mmse_demod = reshape(mmse_demod_,[Nt, 1]);


        cnt_ml = cnt_ml + sum(x_(:,idx)~=ml_demod,"all");
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
title('4 x 4 MIMO, QPSK');


function hat = ml_detector(h, y, Nt)
    qam_table = [-1-1*1j, 1+1*1j, -1+1*1j, 1-1*1j];
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
                    result(cnt) = norm(y - h*xx).^2;
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
    ipHat(find(y_re < 0 & y_im <= 0)) = -1-1*1j;
    ipHat(find(y_re >= 0 & y_im > 0)) = 1+1*1j;
    ipHat(find(y_re < 0 & y_im >= 0)) = -1+1*1j;
    ipHat(find(y_re >= 0 & y_im < 0)) = 1-1*1j;
end