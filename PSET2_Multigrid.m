%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Problem set 2
%Alejandra Torres León
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc

%%%interp(V,nk,gridnew)

%% Parameters
beta = 0.987;
mu = 2;
alpha = 1/3;
delta = 0.012;
rho = 0.95;
sigma = 0.007;
size1 = 100;

%% Value function iteration

%TFP shock
[Tran, sout]=tauchen(rho, sigma, 0, 7, size1, 3); %Use Tauchen function

%TFP grid
nz = exp(sout);

%Capital grid
kss = ((beta*alpha*1)/(1-beta*(1-delta)))^(1/(1-alpha));

%Initial grid
nk= linspace(0.75*kss,1.25*kss,size1)';   %size1 is already predetermined as 100
%% Multigrid


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
   k = nk(argmax(:,:));
norm = max((abs(V(:)-V0(:))));
   if norm>toler
       V0=V;
   else
       break
   end
end
for iz=1:length(nz)
   nk_mg1(:,iz) = interp1(nk,V(:,iz), linspace(0.75*kss,1.25*kss,500)');

end
toc;

%%
%Value function plot
plot(k,V)
xlabel("k'")
ylabel("Value function")
title("Multigrid value function iteration - 100 points")
multigrid1=gcf
saveas(multigrid1, 'multigrid1.png')
%%
%Policy function plot
plot(k)
xlabel("k'")
ylabel("Policy function")
title("Policy function -  100 points")
policy_mg3=gcf
saveas(policy_mg3, 'policy_mg1.png')
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
title("Euler equation error - 100 points")
euler_31=gcf
saveas(euler_31, 'euler_31.png')




%%

%Maximization 
toler  = 10^-5;
maxiter = 100;
nk = linspace(0.75*kss,1.25*kss,500)';
V0 = nk_mg1;

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
    k = nk(argmax(:,:));
norm = max((abs(V(:)-V0(:))));
   if norm>toler
       V0=V;
   else
       break
   end
end
for iz=1:length(nz)
   nk_mg2(:,iz) = interp1(nk,V(:,iz), linspace(0.75*kss,1.25*kss,5000));

end


toc;
%%
%Value function plot
plot(k,V)
xlabel("k'")
ylabel("Value function")
title("Multigrid value function iteration - 500 points")
multigrid2=gcf
saveas(multigrid2, 'multigrid2.png')
%%
%Policy function plot
plot(k)
xlabel("k'")
ylabel("Policy function")
title("Policy function -  500 points")
policy_mg2=gcf
saveas(policy_mg2, 'policy_mg2.png')
%%

% Euler equation error
%Define optimal consumption using k'

c1 = zeros(length(nk),length(nz));
c2 = zeros(length(nk),length(nz));

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
title("Euler equation error - 500 points")
euler_32=gcf
saveas(euler_32, 'euler_32.png')


%%

nk = linspace(0.75*kss,1.25*kss,5000)';
%Maximization 
toler  = 10^-5;
maxiter = 10;

V0 = nk_mg2;

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
    k = nk(argmax(:,:));
norm = max((abs(V(:)-V0(:))));
   if norm>toler
       V0=V;
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
title("Multigrid value function iteration - 5000 points")
multigrid3=gcf
saveas(multigrid3, 'multigrid3.png')

%%
%Policy function plot
plot(k)
xlabel("k'")
ylabel("Policy function")
title("Policy function -  5000 points")
policy_mg3=gcf
saveas(policy_mg3, 'policy_mg3.png')
%%

% Euler equation error
%Define optimal consumption using k'
c1 = zeros(length(nk),length(nz));
c2 = zeros(length(nk),length(nz));

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
title("Euler equation error - 5000 points")
euler_33=gcf
saveas(euler_33, 'euler_33.png')



