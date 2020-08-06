function u = AnnFilter(Z, K)
%% 输入参数：
    %Z              ——测量数据
    %K              ——脉冲个数
%% 输出参数：
    %tk             ——时延参数
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(Z);        % 测量长度
c = Z(1:N-K); 
r = Z(N-K:N).';
M = hankel(c,r);
ZM = M(:,2:end);
y = -1*M(:,1);
% x = (ZM.'*ZM)\ZM.'*y;
x = ZM\y;
A = [1;x];
u = roots(A).^(-1);

end

