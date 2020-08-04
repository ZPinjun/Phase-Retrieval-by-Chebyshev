function [i, j] = UpdateSupp(A, y, x, S, J1, J2)
% 计算待交换的的索引值
% 确定i
All_i = setxor(S,J1);
i1 = find(abs(x(All_i)-min(x(All_i)))<=1e-10);
i = All_i(i1(1));

% 确定j (这里按照一般情况，不利用FFT算法)
All_j = setxor(S,J2);
b = (abs(conj(A.')*x)).^2-y;
Deta_f = abs(A*diag(b,0)*conj(A.')*x);
j1 =  find(abs(Deta_f(All_j)-max(Deta_f(All_j)))<=1e-10);
index = randperm(length(j1));
j = All_j(j1(index(1)));
end


