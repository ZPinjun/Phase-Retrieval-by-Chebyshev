m0 = [0, 4, 6, 12, 20, 42, 50, 59, 71, 100];
mcomb = nchoosek(m0,2);
dm = mcomb(:,2)-mcomb(:,1);
dm(dm>50) = 100 - dm(dm>50);

m = zeros(12,1);
m(1) = 0;
m(end) = 100;

