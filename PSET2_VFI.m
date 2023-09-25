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
[Tran, sout, probst, xt, st]=tauchen(rho, sigma, 0, 7, size1, 3); %Use Tauchen function
plot(st)
z=zeros(0,size1);


%TFP grid
nz = exp(sout);

%Capital grid
kss = ((beta*alpha*1)/(1-beta*(1-delta)))^(1/(1-alpha));
nk = linspace(0.75*kss,1.25*kss,size1)';   



%%
%Maximization 
toler  = 10^-5;
maxiter = 100;

V0 = zeros(length(nk),length(nz));

tic;

for iter = 1:maxiter
for iz=1:length(nz)
for ik=1:length(nk)
c = nz(iz)*(nk(ik).^alpha) + (1-delta)*nk(ik) - nk;
c(c<0)=0;
u = (c.^(1-mu)-1)/(1-mu);
u(c==0) = -1e12;

[V(ik,iz),argmax(ik,iz)] = max(u+ beta*(V0*Tran(:,iz)));

end 
end
    k = nk(argmax);
norm = max((abs(V(:)-V0(:))));
   if norm>toler
       V0=V
   else
       break
   end
end
toc;
%%
%Value function plot
plot(k,V)
xlabel("k'")
ylabel("Value function")
title("Brute force value function iteration")
brute_force=gcf
saveas(brute_force, 'brute_force.png')
%%
%Policy function plot
plot(k)
xlabel("k'")
ylabel("Policy function")
title("Policy function -  Value function iteration")
policy_bf=gcf
saveas(policy_bf, 'policy_bf.png')

%%

% Euler equation error
%Define optimal consumption using k'
for iz = 1:length(nz)
        c1(:,iz)=nz(iz).*(nk.^alpha) + (1-delta).*nk - k(:,iz);
end

for iz = 1:length(nz)
        c2(:,iz)=nz(iz).*(nk.^alpha) + (1-delta).*nk - k(:,iz);
end


%Iterate over the euler equation
for iz = 1:length(nz)
for ik=1:length(nk)
    if ik<length(nk) && iz<length(nz)
       euler_up(ik,:)= beta*(c2(ik+1,:).^(-mu))*(alpha*nz(iz+1).*k(ik+1).^(alpha-1) + (1-delta));
       euler_down(ik,:) = c1(ik,:).^(-mu);
    else 
       euler_up(ik,:)= beta*(c2(ik,:).^(-mu))*(alpha*nz(iz).*k(ik).^(alpha-1) + (1-delta));
       euler_down(ik,:) = c1(ik,:).^(-mu);
end
end
end
euler = euler_up./euler_down;

error = log10(abs(1-euler));
plot(k,error)

xlabel("k'")
ylabel("Log10 Euler Equation Error")
title("Euler equation error")
euler_1=gcf
saveas(euler_1, 'euler_1.png')

