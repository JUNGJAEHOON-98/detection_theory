clear;
N = 10^4; % number of symbols
Eb_N0_dB = [0:3:30]; % multiple Eb/N0 values
Nt = 4;
Nr = 4;

ip = [(2*(rand(1,N)>0.5)-1) + 1j*(2*(rand(1,N)>0.5)-1)];

x_ = reshape(ip, [Nt, N/Nt]);


for Eb_idx = 1:length(Eb_N0_dB)
    disp(Eb_N0_dB(Eb_idx));
    P = sqrt((10^(Eb_N0_dB(Eb_idx)/10))/Nt);
    x = P/sqrt(2) * x_; % normalization of energy to P
    cnt_zfsic = 0;
    cnt_zf = 0;
    cnt_mmse = 0;
    cnt_ml = 0;
    for idx = 1:N/Nt
    % Channel and noise Noise addition
        h = 1/sqrt(2)*[randn(Nr, Nt) + 1j*randn(Nr, Nt)]; % Rayleigh channel
        n = 1/sqrt(2)*[randn(Nr,1) + 1j*randn(Nr,1)];  % white gaussian noise, 0dB variance
        y = h * x(:,idx) + n;


        w_zf = inv(h'* h) * h';
        w_mmse = inv(h'*h + 1/P^2*eye(Nt))*h';

        ml_demod = ml_detector(h, y, P, Nt);

        zf_demod_ = qam_demod(w_zf * y);
        zf_demod = reshape(zf_demod_,[Nt, 1]);

        zfsic_demod = zf_sic(w_zf, h, y, Nt, P);

        mmse_demod_ = qam_demod(w_mmse * y);
        mmse_demod = reshape(mmse_demod_,[Nt, 1]);


        cnt_zfsic = cnt_zfsic + sum(x(:,idx)~=zfsic_demod,"all");
        cnt_ml = cnt_ml + sum(x(:,idx)~=ml_demod,"all");
        cnt_zf = cnt_zf + sum(x_(:,idx)~=zf_demod,"all");
        cnt_mmse = cnt_mmse + sum(x_(:,idx)~=mmse_demod,"all");

    end
    ser_zfsic(Eb_idx) = cnt_zfsic/N;
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
hold on
semilogy(Eb_N0_dB, ser_zfsic, 'd-','Color','#EDB120','LineWidth',2);
legend('ZF', 'MMSE', 'ML', 'ZF-SIC');
xlabel('SNR[dB]')
ylabel('SER');
ylim([10^-3.5 10^0]);
title('4 x 4 MIMO, QPSK');
grid on


function hat = ml_detector(h, y, P, Nt)
    qam_table = P/sqrt(2) * [-1-1*1j, 1+1*1j, -1+1*1j, 1-1*1j];
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
    ipHat(find(y_re < 0 & y_im < 0)) = (-1-1*1j);
    ipHat(find(y_re > 0 & y_im > 0)) = (1+1*1j);
    ipHat(find(y_re < 0 & y_im > 0)) = (-1+1*1j);
    ipHat(find(y_re > 0 & y_im < 0)) = (1-1*1j);
end


function x_hat = zf_sic(w_zf, h, y, Nt, P)
    w_norm1 = vecnorm(w_zf.');
    [B,I1] = mink(w_norm1, 1,'ComparisonMethod','abs');
    aa = w_zf(I1,:)*y;
    x_hat1 = P/sqrt(2) *qam_demod(aa);

    y1 = y-h(:, I1)*x_hat1;
    h(:, I1) = 0;
    w_zf2 = pinv(h'*h)*h';
    w_norm2 = vecnorm(w_zf2.');
    [B,I2] = mink(w_norm2, 2,'ComparisonMethod','abs');
    bb = w_zf2(I2(2),:)*y1;
    x_hat2 = P/sqrt(2) *qam_demod(bb);

    y2 = y1-h(:, I2(2))*x_hat2;
    h(:, I2(2)) = 0;
    w_zf3 = pinv(h'*h)*h';
    w_norm3 = vecnorm(w_zf3.');
    [B,I3] = mink(w_norm3, 3,'ComparisonMethod','abs');
    cc = w_zf3(I3(3),:)*y2;
    x_hat3 = P/sqrt(2) *qam_demod(cc);

    y3 = y2-h(:, I3(3))*x_hat3;
    h(:, I3(3)) = 0;
    w_zf4 = pinv(h'*h)*h';
    w_norm4 = vecnorm(w_zf4.');
    [B,I4] = mink(w_norm4, 4,'ComparisonMethod','abs');
    dd = w_zf4(I4(4),:)*y3;
    x_hat4 = P/sqrt(2) *qam_demod(dd);

    x_hat = zeros([Nt,1]);
    x_hat(I1) = x_hat1;
    x_hat(I2(2)) = x_hat2;
    x_hat(I3(3)) = x_hat3;
    x_hat(I4(4)) = x_hat4;

end