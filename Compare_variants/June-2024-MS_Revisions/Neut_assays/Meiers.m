close all
clear all

lambda = 0.9;
T50 = 4.2;
tau = 65;

Tik = 6.6; % vac -> Delta
Tik = 7.1; % booster ->BA.1
%Tik = 2.2; % vac ->BA.1


delta_t = linspace(0,600);
H_s = zeros(1,length(delta_t));
for i = 1:length(delta_t)
    T = Tik -delta_t(i)./tau;
    H(i) = CalculateH(T,T50,lambda);
end

plot(delta_t,H,'-o')
title('Vac \rightarrow BA.1')
xlabel('Time after exposure \Delta t')
ylabel('c_i^k(\Delta t)')
ylim([0 1])
set(gca,'FontSize',16)

function H = CalculateH(T,T50,lambda)
    H = 1./(1 + exp(-lambda.*(T-T50)));
end