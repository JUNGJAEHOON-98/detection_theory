clear;
P=1;
Nt=4;
qam_table = [-3 -1 1 3];
xx = zeros(Nt, 1);
x_ = zeros(Nt, 1, length(qam_table).^Nt);

T_inv = [-3 0 -8 0; 1 0 3 0; 1 3 0 0; 0 0 1 3];

cnt = 1;
for idx = 1:length(qam_table)
    xx(1) = qam_table(idx);
    for jdx = 1:length(qam_table)
        xx(2) = qam_table(jdx);
        for kdx = 1:length(qam_table)
            xx(3) = qam_table(kdx);
            for hdx = 1:length(qam_table)
                xx(4) = qam_table(hdx);
                aa = T_inv * xx;
                x_1(cnt) = aa(1);
                x_2(cnt) = aa(2);
                x_3(cnt) = aa(3);
                x_4(cnt) = aa(4);
                cnt = cnt + 1;
            end
        end
    end
end
x_1 = unique(x_1);
x_2 = unique(x_2);
x_3 = unique(x_3);
x_4 = unique(x_4);

t = 11.1;
x_1 - t
abs(x_1 - t)
[M I] = min(abs(x_1 - t))
t_hat = x_1(I)