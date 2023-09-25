function [P_t,sout,probst,xt,st]=tauchen(rho,sigma,mu,N,T,m)


if N==1
    P_t=1;
    arho=1;
    sout=0;
    probst=1;
    asigmay=0;

else

%Grid
ub_t = m*sigma/(sqrt(1-rho^2));
lb_t = -ub_t;

grid_t = linspace(lb_t,ub_t,N)


%Transition matrix
P_t = zeros(N,N); %Define an empty matrix of 9x9
delta = (ub_t-lb_t)/(N-1); %Define the half-distance to both sides for a point $z_j$ in the grid

%Loop over the different points in the grid to obtain the associated value
%of the cumulative distribution function
for i=1:length(grid_t)
for j=1:length(grid_t)
    z_i = grid_t(i)
    z_j = grid_t(j)

P_t(i,j)=normcdf((z_j-mu*(1-rho)-rho*z_i+delta/2)/sigma)-normcdf((z_j-mu*(1-rho)-rho*z_i-delta/2)/sigma)

end
end

%Corner points
P_t(1,1)=normcdf((lb_t-mu*(1-rho)-rho*grid_t(1)+delta/2)/sigma) 
P_t(N,N) = 1-normcdf((ub_t-mu*(1-rho)-rho*grid_t(N)-delta/2)/sigma) 

%Calculate discretized time series

    sout=grid_t;

    probst=ones(1,N)/N;
    while 1
        probstnew=probst*P_t;
        lixo=max(abs(probstnew-probst));
        if lixo<1e-9
            break
        end
        probst=probstnew;
    end

end

if ~exist('T','var')
    return
end

rand('seed',0);
epsilon=rand(1,T);
xt=ones(1,T);

xt(1)=floor(N/2)+1;

aux=cumsum(P_t')';
for t=2:T
    [lixo,xt(t)]=max(aux(xt(t-1),:)>epsilon(t));
end
st=grid_t(xt);
    
   

