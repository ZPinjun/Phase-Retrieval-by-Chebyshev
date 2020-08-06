clc;clear;close all
M = 100;
N = 300;
K = 3;
%% --------------------------------------------------- �����ź�
h = [1;5;3];
x = zeros(M,1);
r = 1;
x(r:r+2,1) = 1.2*h;
r = 40;
x(r:r+2,1) = 1.6*h;
r = 69;
x(r:r+2,1) = 2.8*h;

c00 = [1.2; 1.6; 2.8];
m00 = [1; 40; 69];
% m0 = [39;68;29];
y = abs(fft(x,N)).^2;



% %% --------------------------------------------------- ���ݴ���
% H = fft(h,N);
% H = reshape(H,[N,1]);
% xi = y./(abs(H)).^2;

%% --------------------------------------------------- ��֤ƽ���͹�ϵ
% ccomb = nchoosek(c,2);
% mcomb = nchoosek(m,2);
% m_m = mcomb(:,2) - mcomb(:,1);
% 
% left = xi - norm(c,2)^2.*ones(size(xi));
% 
% u_1 = exp(1i*(2*pi/N).*m_m);
% uu = [u_1; conj(u_1)];
% cc = [ccomb(:,1).*ccomb(:,2);ccomb(:,1).*ccomb(:,2)];
% 
% right = zeros(N,1);
% for i = 1:N
%     right(i) = sum(cc.*(uu.^(i-1)));
% end
% left - right % ������ȣ�˵���Ƶ�û�д���

%% --------------------------------------------------- �����㻯�˲�����
% tkrec = AnnFilter(left, K*(K-1), 1).';  % ����̫��

% %% --------------------------------------------------- �Ľ����б�ѩ��
% % ���ݴ���
% H = fft(h,N);
% H = reshape(H,[N,1]);
% xi = y./(abs(H)).^2;
% %  ��ʼ��
% threshold = 10e-4;
% T = 50;
% beta = 1;
% 
% ccomb = nchoosek(c0,2);
% c1 = 2*ccomb(:,1).*ccomb(:,2);
% % m1 = [39;68;29]; m1 = [42;70;30];
% % m1 = ones(3,1)*(M/3);
% m1 = floor(rand(floor(K*(K-1)*0.5),1)*M); 
% f1_m = @(mx) (-2*pi/N)*diag(0:N-1,0)*sin((2*pi/N)*(0:N-1).'*mx.')*diag(c1,0);
% f_m = @(mx) cos((2*pi/N)*(0:N-1).'*mx.')*c1-xi+ones(N,1)*norm(c0,2)^2;
% % ����
% err = zeros(T+1,1);
% alpha = 0.001;
% disp(['��',num2str(1),'�ε�����~',', m=',num2str(m1.')]);
% for tt = 1:T
% % ------------------------ �ݶ��½���
%     dd = alpha*2*f1_m(m1).'*f_m(m1);
%     m1 = m1 - dd;
%     err(tt+1) = (2*K)\norm(f_m(m1),2)^2;
%     disp(['��',num2str(tt+1),'�ε�����',num2str(err(tt+1)),', m=',num2str(m1.'), ', ����',num2str(dd.')]);
% % ------------------------ �б�ѩ�򲻺�ʹ    
% %     disp(['��',num2str(tt),'�ε�����~',', m=',num2str(m1.')]);
% %     f1m = f1_m(m1);
% %     gm = (f1m.'*f1m)\f1m.';
% %     fm = f_m(m1);
% %     y = m1 - beta*gm*fm;
% %     dd = gm*(fm+beta^2*(f_m(y)-fm+beta*f1m*gm*fm));
% %     m1 = m1 - gm*(fm+beta^2*(f_m(y)-fm+beta*f1m*gm*fm));
% % %     m1 = floor(m1);
% %     err(tt+1) = (2*K)\norm(fm,2)^2;
% %     disp(['��',num2str(tt+1),'�ε�����',num2str(err(tt+1)),', m=',num2str(m1.'), ', ����',num2str((gm*(fm+beta^2*(f_m(y)-fm+beta*f1m*gm*fm))).')]);
% %     figure(10)
% %     plot(1:tt,err(2:tt+1));
% end
% 
% m = round(m1);
% err(tt+2) = (2*K)\norm(f_m(m),2)^2;
% disp(['��',num2str(tt+2),'�ε�����',num2str(err(tt+2))]);



%% --------------------------------------------------- get xi_hat
H = fft(h,N);
H = reshape(H,[N,1]);
xi = y./(abs(H)).^2-ones(N,1)*norm(c00,2)^2;

mcomb = nchoosek(m00,2);
m1 = (2*pi/N)*(mcomb(:,2)-mcomb(:,1));
ccomb = nchoosek(c00,2);
c1 = 2*ccomb(:,1).*ccomb(:,2);
left = zeros(N,1);
for i = 1:N
    left(i,1) = c1.'*(cos(m1).^(i-1));
end

% ���㲢��������n����ϵ����
% maxn = 1000;
% cosn_coeff = zeros(maxn,round(maxn/2));
% for n = 1:maxn
%     for k = 0:floor(n/2)
%         if(k~=0)
%             cosn_coeff(n,k+1) = 2^(-2*k)*(nchoosek(n-k,k)+nchoosek(n-k-1,k-1))*(-1)^(k+1);
%         else
%             cosn_coeff(n,k+1) = 2^(1-n);    
%         end
%     end
% end
% save('cosn_coeff.mat','cosn_coeff');
load('cosn_coeff.mat');
% calculate xi_hat
for j = 3:N
    if mod(j,2)==1
        xi1 = xi(j:-2:mod(j,2));
    else
        xi1 = xi(j:-2:mod(j,2)+2);
    end
    xi(j) = cosn_coeff(j-1,1:floor((j+1)/2))*xi1;
end

% call prony method
u = AnnFilter(xi(1:100), 3);

u
cos(m1)
dmr = acos(u).*(N/(2*pi))