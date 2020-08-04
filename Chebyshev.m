function [cr, err] = Chebyshev( y, h, m )
% Chebyshev法迭代优化幅度参数
%% --------------------------------------------------- 数据处理
K = length(m);
N = length(y);
H = fft(h,N);
H = reshape(H,[N,1]);
xi = y./(abs(H)).^2;
xi(1) = sqrt(xi(1));

%% --------------------------------------------------- 初始化
threshold = 10e-4;
T = 12;
beta = 1;
es = (sum(y))/(K*N*norm(h,2)^2);
c0 = es*ones(K,1);
% load('c0.mat');

n = (1:N-1).';
mcomb = nchoosek(m,2);
q = mcomb(:,1) - mcomb(:,2);
term1 = 2*cos((2*pi/N)*n*q.');
p = @(a) (nchoosek(a,2)*[1;0]).*(nchoosek(a,2)*[0;1]);
f = @(a) [sum(a); norm(a,2)^2*ones(N-1,1)+term1*p(a)]-xi;

s1 = kron(ones(N-1,1),m)*ones(1,K) - ones((N-1)*K,1)*m.';
s2 = kron(n,ones(K,1))*ones(1,K);
term2 = cos((2*pi/N)*(s1).*(s2));

err = zeros(T,1);
c = c0;
%% --------------------------------------------------- 迭代
for i = 0:T-1
    fc = f(c);
    f1c = [ones(1,K);2*kron(eye(N-1),c.')*term2];
    gc = (f1c.'*f1c)\f1c.';
    z = c - beta*gc*fc;
%     disp(['第',num2str(i+1),'次幅度：',num2str(c.')]);
%     disp(['第',num2str(i+1),'次方向：',num2str((gc*(fc+beta^2*(f(z)-fc + beta*f1c*gc*fc   ))).')]);
    c = c - gc*(fc+beta^2*(f(z)-fc + beta*f1c*gc*fc   ));
    err(i+1) = (2*K)\norm(f(c),2)^2;
    disp(['第',num2str(i+1),'次迭代误差：',num2str(err(i+1))]);
    figure(10)
    plot(0:i,err(1:i+1));
end
 
cr = c;

end

