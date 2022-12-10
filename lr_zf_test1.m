clear;

P = 1; % 일단 테스트용이라 임의로 설정
qam_table = [-1-1*1j, 1+1*1j, -1+1*1j, 1-1*1j];
N = 2;
Nt = 2;
Nr = 2;

h = 1/sqrt(2)*[randn(Nr, Nt) + 1j*randn(Nr, Nt)];
n = 1/sqrt(2)*[randn(Nr,1) + 1j*randn(Nr,1)];
ip = [(2*(rand(1,N)>0.5)-1) + 1j*(2*(rand(1,N)>0.5)-1)];
x = reshape(ip, [Nt, N/Nt]);

xx = zeros(Nt, 1);

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
                aa = pinv(T) * xx;
                z1_set_(cnt) = aa(1);
                z2_set_(cnt) = aa(2);
                cnt = cnt + 1;
    end
end
z1_set_
z1_set = unique(z1_set_)
z2_set = unique(z2_set_);


[M I] = min(abs(z_hat(1) - z1_set));
z_hat(1) = z1_set(I);
[M I] = min(abs(z_hat(2) - z2_set));
z_hat(2) = z2_set(I);


x_hat = T*z_hat;



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