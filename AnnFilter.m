function tk = AnnFilter(x, K, T)
%% OMP重构算法_仿真实验_现代采样理论、技术及应用
%% 输入参数：
    %x              ——傅里叶系数
    %K              ——脉冲个数
    %T              ——周期数
%% 输出参数：
    %tk             ——脉冲时延
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(x);        % 获取傅里叶系数个数
M = floor(N/2);       % 傅里叶系数个数的一半 
L = K;                % 滤波器长度

r = x( (-M+L: -1: -M) + floor(N/2)+1 );   % 将傅里叶系数截取后进行重新排列
c = x( (-M+L: M) + floor(N/2) + 1 );     % 将傅里叶系数截取后进行重新排列
my_A = toeplitz(c, r);   % 利用重新排列后的傅里叶系数矩阵构造Topelitz矩阵

[~, ~, V] = svd(my_A, 0);    % 对构造的Topelitz矩阵进行奇异值分解

ann_filter_coeff = V(:, K+1);% 选择最小的奇异值所对应的右奇异向量作为滤波器系数向量

tk = sort(mod(-angle(roots(ann_filter_coeff))/(2*pi)*T, T)); % 通过最小二乘法恢复脉冲时延
