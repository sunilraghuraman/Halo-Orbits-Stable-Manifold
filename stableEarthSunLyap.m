clear all
clc
format long;
epsilon=2*10^(-8);
% haloPeriod= 2.69239959528586;
haloPeriod= 3.0562;
T=3.1*haloPeriod;
k=20;
G=1;
% GM_earth=3.986*10^(5);
% GM_sun=1.327*10^(11);   
au=1.496*10^8;
% mu=GM_earth/(GM_moon+GM_earth);
mu = 0.0121505856;
m1=1-mu;
m2=mu;
r0 = [1-mu; 0; 0]+(1/au)*[-1200000; 0; -280000];
v0 = (1/au)*(60*60*24*365.25)/(2*pi)*[0; -0.350; 0];
C=jacobiConst(r0,v0,mu);
Lagrange_Points


% outL1=L1;
% outL2=L2;
% outL3=L3;
% outL4=L4;
% outL5=L5;
% v_equil=[0; 0; 0];
% 
% %Compute the energy at each libration point
% CL3=jacobiConst(L3, v_equil, mu);
% CL1=jacobiConst(L1, v_equil, mu);
% CL2=jacobiConst(L2, v_equil, mu);
% CL4=jacobiConst(L4, v_equil, mu);
% CL5=CL4; 

t0=0;
numSteps=300;
tspan=linspace(0,haloPeriod+0.1*haloPeriod,numSteps);
y0=[r0; v0]';                                  
options=odeset('RelTol',1e-13,'AbsTol',1e-22);  
[t,x_halo]=ode45('CRTBP',tspan,y0,options,[],G,mu);

%compute the monodromy matrix
Phi=stateTransCRTBP(t0, haloPeriod, [r0; v0], mu);

%compute the eigenvalues and eigenvectors

%Get eigenvalues and eigenvectors
[Phi_eigVecs, Phi_eigValOnDiag]=eigs(Phi);
PhiEigs=diag(Phi_eigValOnDiag);

figure(1)
hold on
h=round((numSteps-1)/k);
for n=1:k
haloPoint(n,1:6)=x_halo((n-1)*h+1,1:6);
haloTimes(n)=t((n-1)*h+1);
end




I=eye(6);
%The un-stable eigenvector is
deltaU=Phi_eigVecs(1:6,6);

%Compute an orbit (fiber of the unstable manifold) on the unstable 
%manifold begining at the origin of the orbit.  first the initial
%condition is found:
xU_p(1,1:6)=haloPoint(1,1:6)'+epsilon*I*deltaU;


%compute the stable orbit 
t0=0;
numSteps=2000;
tspan=linspace(0, T,numSteps);
%numerical integration
y0=xU_p(1,:)';                                      %inital condition
options=odeset('RelTol',2.5e-13,'AbsTol',1e-22);     %set tolerences
[t,xU_p1]=ode45('BackwardCRTBP',tspan,y0,options,[],G,mu);




for n=2:k
    %compute the initial conditions for the manifold orbits.  These are
    %found by pushing the vector 'deltaU' around with the state
    %transition matrix.  (This is a push forward of deltaU). The
    %justification of this is in the 5th note set.
    initialconditions=n;
    perturbationVector=stateTransCRTBP(t0, haloTimes(n), [r0; v0], mu)*deltaU;
    u_pV=perturbationVector/norm(perturbationVector);
    xU_p(n,1:6)=haloPoint(n,1:6)'+epsilon*u_pV;
end



%Compute and plot the unstable manifold

for n=2:k
    manifoldLoop=n;
    t0=0;
    numSteps=2000;
    tspan=linspace(0, T,numSteps);
    y0=xU_p(n,1:6)';                                     
    options=odeset('RelTol',2.5e-13,'AbsTol',1e-22);     
    [t,xU_pOrb]=ode113('BackwardCRTBP',tspan,y0,options,[],G,mu);

    %plot the unstable pos orbits as red lines
    plot3(xU_pOrb(:,1), xU_pOrb(:,2), xU_pOrb(:,3), 'r')
end

%plot the libration points
% plot3(L1(1), L1(2), L1(3), 'k*', L2(1),  L2(2), L2(3), 'k*')


%plot the planets/primaries 
plot3(1-mu, 0, 0,'g*')
plot3(-mu, 0, 0, 'b*')

%plot the periodic orbit in black dots
plot3(x_halo(:,1), x_halo(:,2), x_halo(:,3), 'k.')
