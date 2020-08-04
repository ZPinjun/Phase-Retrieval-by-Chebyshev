function [i, j] = UpdateSupp(A, y, x, S, J1, J2)
% ����������ĵ�����ֵ
% ȷ��i
All_i = setxor(S,J1);
i1 = find(abs(x(All_i)-min(x(All_i)))<=1e-10);
i = All_i(i1(1));

% ȷ��j (���ﰴ��һ�������������FFT�㷨)
All_j = setxor(S,J2);
b = (abs(conj(A.')*x)).^2-y;
Deta_f = abs(A*diag(b,0)*conj(A.')*x);
j1 =  find(abs(Deta_f(All_j)-max(Deta_f(All_j)))<=1e-10);
index = randperm(length(j1));
j = All_j(j1(index(1)));
end


