clear all;
close all;
N = 10^5; % number of symbols
Eb_N0_dB = [0:3:30]; % multiple Eb/N0 values
Nt = 4;
Nr = 4;

ip = [(2*(rand(1,N)>0.5)-1) + 1j*(2*(rand(1,N)>0.5)-1)];
x_ = reshape(ip, [Nt, N/Nt]);

for Eb_idx = 1:length(Eb_N0_dB)
    disp(Eb_N0_dB(Eb_idx));
    P = sqrt((10^(Eb_N0_dB(Eb_idx)/10))/Nt);
    x = P/sqrt(2) * x_; % normalization of energy to P

    cnt_sd = 0;
    cnt_sdr = 0;
    cnt_lrzf = 0;
    cnt_mmsesic = 0;
    cnt_zfsic = 0;
    cnt_zf = 0;
    cnt_mmse = 0;
    cnt_ml = 0;
    tic()
    for idx = 1:N/Nt
    % Channel and noise Noise addition
        h = 1/sqrt(2)*[randn(Nr, Nt) + 1j*randn(Nr, Nt)]; % Rayleigh channel
        n = 1/sqrt(2)*(randn(Nr,1) + 1j*randn(Nr,1));  % white gaussian noise, 0dB variance
        y = h * x(:,idx) + n;


        w_zf = inv(h'* h) * h';
        w_mmse = inv(h'*h + 1/P^2*eye(Nt))*h';

        sd_demode = sd_detector(h, x(:,idx), n, Nt, P);
        sdr_demode = sdr_detector(h, x(:,idx) ,n, P, Nt);
        lrzf_demod = lr_zf(x_(:,idx), h, n, P, Nt);

        ml_demod = ml_detector(h, y, P, Nt);

        zf_demod_ = qam_demod(w_zf * y);
        zf_demod = reshape(zf_demod_,[Nt, 1]);

        mmsesic_demod = mmse_sic(w_mmse, h, y, Nt, P);
        zfsic_demod = zf_sic(w_zf, h, y, Nt, P);

        mmse_demod_ = qam_demod(w_mmse * y);
        mmse_demod = reshape(mmse_demod_,[Nt, 1]);


        cnt_sd = cnt_sd + sum(x_(:,idx)~=sd_demode,"all");
        cnt_sdr = cnt_sdr + sum(x(:,idx)~=sdr_demode,"all");
        cnt_lrzf = cnt_lrzf + sum(round(x(:,idx), 10)~=round(lrzf_demod, 10),"all");
        cnt_mmsesic = cnt_mmsesic + sum(x(:,idx)~=mmsesic_demod,"all");
        cnt_zfsic = cnt_zfsic + sum(x(:,idx)~=zfsic_demod,"all");
        cnt_ml = cnt_ml + sum(x(:,idx)~=ml_demod,"all");
        cnt_zf = cnt_zf + sum(x_(:,idx)~=zf_demod,"all");
        cnt_mmse = cnt_mmse + sum(x_(:,idx)~=mmse_demod,"all");

    end
    toc()
    ser_sd(Eb_idx) = cnt_sd / N;
    ser_sdr(Eb_idx) = cnt_sdr / N;
    ser_lrzf(Eb_idx) = cnt_lrzf / N;
    ser_mmsesic(Eb_idx) = cnt_mmsesic/N;
    ser_zfsic(Eb_idx) = cnt_zfsic/N;
    ser_ml(Eb_idx) = cnt_ml/N;
    ser_zf(Eb_idx) = cnt_zf/N;
    ser_mmse(Eb_idx) = cnt_mmse/N;
end

figure
semilogy(Eb_N0_dB, ser_ml, '-','Color','#000000','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_sd, 'o-','Color', '#FF0000','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_sdr, 'x-','Color', '#0000FF','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_lrzf, '^-','Color','#77AC30','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_zfsic, 'd-','Color','#EDB120','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_mmsesic, 'x-','Color','#4DBEEE','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_zf, 'v--','Color','#EDB120','LineWidth',2);
hold on
semilogy(Eb_N0_dB, ser_mmse, '^--','Color','#4DBEEE','LineWidth',2);

legend('ML', 'SD with d^2=9', 'SDR', 'LR-ZF', 'ZF-SIC', 'MMSE-SIC', 'ZF', 'MMSE');
xlabel('SNR[dB]')
ylabel('SER');
ylim([10^-3.5 10^0]);
title('4 x 4 MIMO, QPSK');
grid on

function result = sd_detector(h, x, n, Nt, P)
    y = (h * x + n)/ (P/sqrt(2));
    h_re_im = [real(h) -imag(h); imag(h) real(h)];
    y_re_im = [real(y); imag(y)];
    xx = zeros(2*Nt,1);

    w_zf = pinv(h_re_im);

    G = h_re_im.'*h_re_im;
    R = chol(G);
    x_hat = w_zf*y_re_im;

    result_x = x_hat;

    up_bound = zeros(1,2*Nt);
    d1_square = zeros(1,2*Nt); % d_k
     % s_k

    x_hat_k1_k2 = zeros(1,2*Nt); % s_^_k|k+1
    caseno = 1;

    while (caseno~=0)
        switch (caseno)
            case 1
                k = 2*Nt;
                d_square = 9;
                d1_square(:,k) = d_square;
                max_distance = -Inf;
                x_hat_k1_k2(k) = x_hat(Nt*2);
                caseno = 2;
            case 2
                if d1_square(:,k)<0
                    up_bound(:,k) = 1;
                    low_bound =1;
                else
                    up_bound(:,k) = (sqrt(d1_square(:,k))/R(k,k)) + x_hat_k1_k2(:,k);
                    low_bound = -sqrt(d1_square(:,k))/R(k,k) + x_hat_k1_k2(:,k);
                    up_bound(:,k) = 1-2*(up_bound(:,k)<1);
                    low_bound = 1-2*(low_bound<-1)-2;
                end
                xx(k,:) = low_bound;
                caseno = 3;
            case 3
                xx(k,:) = xx(k,:) + 2;
                if xx(k,:) <= up_bound(:,k)
                    caseno = 5;
                else
                    caseno = 4;
                end
            case 4
                k = k + 1;
                if k == 2*Nt +1
                    break;
                else
                    caseno = 3;
                end
            case 5
                if k == 1
                    caseno = 6;
                else
                    k = k - 1;
                    d1_square(:,k) = d1_square(:,k+1) - (R(k+1,k+1) * (xx(k+1,:) - x_hat_k1_k2(:,k+1)))^2;
                    x_hat_k1_k2(:,k) = x_hat(k) - (1/R(k,k)) * ( R(k,k+1:end) * ( xx(k+1:end,:) - x_hat(k+1:end) ) );
                    caseno = 2;
                end
            case 6
                tmp_distance = d1_square(:,1)-(R(k,k)*(xx(k)-x_hat_k1_k2(k)))^2;

                if tmp_distance >= max_distance
                    max_distance = tmp_distance;
                    result_x = xx;
                end
                caseno = 3;
        end
    end
    %test_1 = (x_hat>0)-(x_hat<0);
    %R_inv_distance = 9 - sum((R*(test_1 -x_hat)).^2);
    %result_distance = 9 - sum((R*(result_x-x_hat)).^2);
    real_val = (result_x(1:Nt)>0)-(result_x(1:Nt)<0);
    imag_val = (result_x(Nt+1:2*Nt)>0)-(result_x(Nt+1:2*Nt)<0);
    result = real_val +1j*imag_val;

end

function x_hat = sdr_detector(h, x ,n, P, Nt)
    y = (h * x + n)/ (P/sqrt(2));
    y_re_im = [real(y);imag(y)];
    h_re_im = [real(h) -imag(h); imag(h) real(h)];
    L = [h_re_im.'*h_re_im -h_re_im.'*y_re_im; -y_re_im.'*h_re_im y_re_im.'*y_re_im];

    cvx_begin quiet;
        variable X(2*Nt+1, 2*Nt+1);
        minimize(trace(L*X));
        diag(X)==1;
        X == semidefinite(2*Nt+1);
    cvx_end;

    [V, D] = eig(X);
    [M I] = max(diag(D));
    a = V(:, I);
    if a(length(V)) < 0
       a = -a;
    end
    demode_ = a(1:Nt)+1j*a(Nt+1:2*Nt);
    x_hat = reshape((P/sqrt(2)) *qam_demod(demode_), [Nt, 1]);

end

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

function x_hat = mmse_sic(w_mmse, h, y, Nt, P)
    x_hat = zeros(Nt,1);

    for i = 1:Nt
        w_mmse = pinv(h'*h + 1/P^2*eye(Nt))*h';

        signal_power = vecnorm(diag(w_mmse*h), 2, 2).^(2);
        interference = vecnorm(w_mmse*h - diag(diag(w_mmse*h)), 2, 2).^(2);
        w_norm = vecnorm(w_mmse, 2, 2);

        sinr = (signal_power.*(P/sqrt(2))) ./ (interference.*(P/sqrt(2)) + w_norm);
        [B,I] = max(sinr);

        x_hat1 = P/sqrt(2) * qam_demod(w_mmse(I,:) * y);
        x_hat(I) = x_hat1;
        y = y - h(:, I) * x_hat1;
        h(:, I) = 0;
    end
end

function x_hat = lr_zf(x, h, n, P, Nt)
    qam_table = P/sqrt(2) * [-1-1*1j, 1+1*1j, -1+1*1j, 1-1*1j];
    xx = zeros(Nt, 1);
    z1_set_ = zeros(Nt, 1, length(qam_table).^Nt);

    T = MLLL(h, Nt);
    z_ = pinv(T)*x;
    z =  P/sqrt(2) * z_;
    h_til = h*T;

    y = h_til*z + n;

    z_hat = pinv(h_til) * y;

    cnt = 1;

    for idx = 1:length(qam_table)
         xx(1) = qam_table(idx);
         for jdx = 1:length(qam_table)
             xx(2) = qam_table(jdx);
             for kdx = 1:length(qam_table)
                 xx(3) = qam_table(kdx);
                 for hdx = 1:length(qam_table)
                     xx(4) = qam_table(hdx);
                     z = pinv(T) * xx;
                     z1_set_(:,:,cnt) = z;
                     result(cnt) = sum(abs(z_hat - z));
                     cnt = cnt + 1;
                end
            end
        end
    end

    [M I] = min(result);
    z_hat = z1_set_(:, I);

    x_hat = T*z_hat;
end

function T = MLLL(H,Tx)
[Q,R] = qr(H);
%initialize
lambda = 3/4;
mm = Tx;
T     = eye(mm);
S_flag = 0;
iter = 0;
% Set the number of iterations
iter_max = 5;
while S_flag == 0 && iter<iter_max
    S_flag = 1;
    iter = iter+1;
  %size reduction
for k = 2:mm
   for l=k-1:-1:1
      mu = round(R(l,k)/R(l,l));
      if abs(mu) ~= 0
         R(1:l,k) = R(1:l,k) - mu * R(1:l,l);
         T(:,k)   = T(:,k)   - mu * T(:,l);
      end
   end

    temp_div = norm(R(k-1:k,k));
   %Siegel Condition
   if lambda*abs(R(k-1,k-1)) > abs(R(k,k))

      %exchange columns k and k-1
        temp_r=R(:,k-1);
        temp_t=T(:,k-1);

        R(:,k-1)=R(:,k);
        R(:,k)=temp_r;
        T(:,k-1)=T(:,k);
        T(:,k)=temp_t;
      % Given's rotation
      alpha     = R(k-1,k-1) / temp_div;
      beta     = R(k,k-1)   / temp_div;
      Theta = [alpha' beta; -beta alpha];

      R(k-1:k,k-1:end) = Theta * R(k-1:k,k-1:end);
      S_flag = 0;
   end
      k = k+1;
end
end
end