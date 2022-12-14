
% MIMO sphere decoder
% works for 2-PAM only (it can be very easily modified for higher PAM)

clear all
close all
clc
Nt = 4; % number of transmit antennas
Nr = 4; % number of receive antennas
fade_var = 1; % fade variance of the MIMO channel
SNR_dB = 10; % SNR per bit (dB)
epsilon = 0.000001; % for choosing radius
num_frames = 10^2; % simulation runs

% SNR parameters
P = 10^(0.1*SNR_dB)/Nt;
%noise_var = 1*fade_var*Nt*Nr/(2*SNR*1*Nt); %awgn variance

% radius (square) of the sphere
d1_2 = 9;

C_Ber=0;
tic()
for i1=1:num_frames
% source
a = randi([0 1],1,Nt);

% 4-QAM
x = P/sqrt(2) * [(2*(rand(1,Nt)>0.5)-1) + 1j*(2*(rand(1,Nt)>0.5)-1)];

% (Nr x Nt) MIMO flat channel
h = 1/sqrt(2)*[randn(Nr, Nt) + 1j*randn(Nr, Nt)];

% AWGN
noise = 1/sqrt(2)*[randn(Nr,1) + 1j*randn(Nr,1)];

% channel output
chan_op = h*x.' + noise;

%--------------- sphere decoder-----------------------------------
[Q,R1] = qr(h); % QR decomposition
seq2 = zeros(1,Nt); % s_k
R = R1(1:Nt,:);% R matrix
Q1 = Q(:,1:Nt);
Q2 = Q(:,Nt+1:end);
y = Q1'*chan_op;
s_hat = inv(R)*Q'*y;
dec_seq = zeros(1,Nt); % decoded sequence
up_bound = zeros(1,Nt); % upper bound
d2 = zeros(1,Nt);
y1 = zeros(1,Nt); % yk|k+1
caseno = 1;

while (caseno~=0)
switch (caseno)
    case 1
          k = Nt;
          d2(k) = d1_2 - norm(chan_op)^2 - norm(h*s_hat)^2;
          % y1(k) = y(Nt);
          s_hat1(k) = s_hat(Nt);
          caseno=2;
    case 2
         temp = [(-sqrt(d2(k))+y1(k))/R(k,k) (sqrt(d2(k))+y1(k))/R(k,k)];
         up_bound(k) = max(temp) ;
         low_bound = min(temp);
         seq2(k) = low_bound-1;
         caseno=3;
    case 3
         seq2(k) = seq2(k)+1;
         if abs(seq2(k))<=abs(up_bound(k))
             caseno = 5;
         else
             caseno = 4;
         end
    case 4
         k = k+1;
         if k==Nt+1
             break;
         else
             caseno = 3;
         end
    case 5
        if k==1
            caseno = 6;
        else
            k=k-1;
            d2(k) = d2(k+1)-(y1(k+1) - R(k+1,k+1)*seq2(k+1))^2;
            y1(k) = y(k)-R(k,k+1:Nt)*seq2(k+1:end).';
            caseno = 2;
        end
    case 6
        dec_seq=seq2; % dec_seq is decoded sequence
        caseno =3;

end % for switch
end % for while

x
qam_demod(dec_seq)

C_Ber = C_Ber +nnz(dec_seq-x);
end
toc()


function ipHat = qam_demod(input)
    y_re = real(input);
    y_im = imag(input);
    ipHat(find(y_re < 0 & y_im < 0)) = (-1-1*1j);
    ipHat(find(y_re > 0 & y_im > 0)) = (1+1*1j);
    ipHat(find(y_re < 0 & y_im > 0)) = (-1+1*1j);
    ipHat(find(y_re > 0 & y_im < 0)) = (1-1*1j);
end
