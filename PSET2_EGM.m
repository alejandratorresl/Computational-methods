%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem set 2
%Alejandra Torres León
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

%% Parameters
beta = 0.987;
mu = 2;
alpha = 1/3;
delta = 0.012;
rho = 0.95;
sigma = 0.007;
size1 = 500;

%% Value function iteration

%TFP shock
[Tran, sout]=tauchen(rho, sigma, 0, 7, size1, 3); %Use Tauchen function

%TFP grid
nz = exp(sout);

%Capital grid
kss = ((beta*alpha*1)/(1-beta*(1-delta)))^(1/(1-alpha));
nk = linspace(0.75*kss,1.25*kss,size1)';   

toler  = 10^-5;
maxiter = 1;

%% FOC
tic;
for iter=1:maxiter
for iz=1:length(nz)
for ik=1:length(nk)

c0(ik,iz)=alpha.*nz(iz).*nk(ik).^alpha+(1-delta).*nk(ik);


b= @(k) (alpha*nz(iz)*k.^alpha+(1-delta)*k-nk(ik) - (1./(alpha*nz(iz)*k.^alpha+(1-delta)*k-nk(ik)*(beta*(c0.^(-mu))*(alpha*nz(iz).*nk(ik).^(alpha-1) + (1-delta))).^mu)))
K(ik,iz) = fsolve(b,1);


end
end
c1=nz(iz).*(nk.^alpha) + (1-delta).*nk - K

norm = max((abs(c1(:)-c0(:))));
   if norm>toler
       c0=c1;
   else
       break
   end
end
toc;

%%
%Policy function plot
plot(K)
xlabel("k'")
ylabel("Policy function")
title("Policy function -  EGM")
congm=gcf
saveas(congm, 'policy_egm.png')

%%
%Policy function plot
plot(K, c1)
xlabel("k'")
ylabel("Updated consumption policy")
title("Consumption policy -  EGM")
conegm=gcf
saveas(conegm, 'conegm.png')
%%

% Euler equation error
%Define optimal consumption using k'
c1 = zeros(length(nk),length(nz));
c2 = zeros(length(nk),length(nz));

for iz = 1:length(nz)
        c1(:,iz)=nz(iz).*(nk.^alpha) + (1-delta).*nk - K(:,iz);

end

for iz = 1:length(nz)
        c2(:,iz)=nz(iz).*(nk.^alpha) + (1-delta).*nk - K(:,iz);

end

%Iterate over the euler equation
for iz = 1:length(nz)
for ik=1:length(nk)
    if ik<length(nk) && iz<length(nz)
       euler_up(ik,:)= beta*(c2(ik+1,:).^(-mu))*(alpha*nz(iz+1).*K(ik+1).^(alpha-1) + (1-delta));
       euler_down(ik,:) = c1(ik,:).^(-mu);
    else 
       euler_up(ik,:)= beta*(c2(ik,:).^(-mu))*(alpha*nz(iz).*K(ik).^(alpha-1) + (1-delta));
       euler_down(ik,:) = c1(ik,:).^(-mu);
end
end
end
euler = euler_up./euler_down;

error = log10(abs(1-euler));
plot(K,error)

xlabel("k'")
ylabel("Log10 Euler Equation Error")
title("Euler equation error - 5000 points")
euler_egm=gcf
saveas(euler_egm, 'euler_egm.png')



