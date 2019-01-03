%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AEROSP 548 
%Final Project
%Team: Sunil Raghuraman, Monica A. Salunkhe
%Code Description: This code finds an orbit with a given initial position close to a halo orbit using Newton's method.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all 
clc
format long;

G = 1;
mu_sun = 1.327e11;
m_u = 0.0121505856;
pos = 35786000;      %Altitude of Earth parking orbit. 
r_0 = [ 0.8234; 0; 0];
v_0 = [ 0; 0.1263; 0];
er_0 = [pos;0;0];
ev_0 = [0;sqrt(mu_sun/pos);0];
x_0 = [r_0; v_0];
e_0 = [er_0;ev_0];
t_initial = 2.7430;
C_energy_jacobi = Jacobi_Const(r_0, v_0, m_u);
t0 = 0;
t_step = 100;
t_final = 365.25*86400;
t_span = t0:t_step:t_final;
options = odeset('RelTol',1e-12,'AbsTol',1e-6);    
[~,Xp] = ode45(@satdyn,t_span,e_0,options);

%Initial Guess
x_n = x_0;
tau_n = t_initial;
tau = tau_n;

a = 14;
for k = 1:a
    k_out = k;  
    t_0 = 0;
    t_f = tau_n;
    n_steps = 2000;
    t_span = linspace(0, t_f, n_steps);
    y_0 = x_n';                                         
    options = odeset('RelTol',2.5e-14,'AbsTol',1e-22);    
    [~,x_ref] = ode45('CRthreeBP',t_span,y_0,options,[],G,m_u);
    Phi = STM_2_CRTBP(t_0, t_f, x_n, m_u);
    f_x = CRTBP_vec_field(x_ref(n_steps,:), G, m_u);
    %Differential Components
    diff_comp_11 = Phi(4,3);
    diff_comp_12 = Phi(4,5);
    diff_comp_21 = Phi(6,3);
    diff_comp_22 = Phi(6,5);
    diff_comp_31 = Phi(2,3);
    diff_comp_32 = Phi(2,5);
    %Finding differential of the constraint vector
    derv_f = [diff_comp_11, diff_comp_12, f_x(4); diff_comp_21, diff_comp_22, f_x(6); 
              diff_comp_31, diff_comp_32, f_x(2)];
    inv_derv_f = inv(derv_f);
    %Next Target Initial Conditions
    T_1 = [x_n(3); x_n(5); tau_n];
    T_2 = -inv_derv_f*[x_ref(n_steps, 4); x_ref(n_steps, 6); x_ref(n_steps, 2)];
    x_st = T_1 + T_2;
    x_n = [r_0(1); 0; x_st(1); 0; x_st(2); 0];
    tau_n = x_st(3);
end    

C_energy_target_orbit = Jacobi_Const(x_n(1:3), x_n(4:6), m_u);
Delta_initial_energy = abs(C_energy_jacobi - C_energy_target_orbit);

t_0 = 0;
halo_period = 2 * tau_n;
n_steps = 2000;
t_span = linspace(0, halo_period, n_steps);
y_0 = x_n';                                             
options = odeset('RelTol',1e-12,'AbsTol',1e-12);    
[~,x_halo]=ode45('CRthreeBP',t_span,y_0,options,[],G,m_u);

periodicity = norm(x_halo(1,1:3) - x_halo(n_steps,1:3)); %Order of magnitude to which the orbit is periodic
Phi = STM_2_CRTBP(t_0, halo_period, x_n, m_u);
[Phi_eigen_vectors, Phi_eigen_values] = eigs(Phi);

t_f = 2 * t_f;
n_steps = 2000;
t_span = linspace(0, t_f, n_steps);
y_0 = x_0';                                           
options = odeset('RelTol',1e-12,'AbsTol',1e-22);   
[t,x] = ode45('CRthreeBP',t_span,y_0,options,[],G,m_u);

h = Lagrange_Points_2;

figure(1)
hold on
plot(Xp(:,1)/384400000, Xp(:,2)/384400000,'k','linewidth',1);
[X,Y] = meshgrid(-1.1:0.002:1.2);
r_1 = ((X + m_u).^2+Y.^2).^(1/2);
r_2 = ((X + m_u-1).^2+Y.^2).^(1/2);
P = 2*((1/2)*(X.^2 + Y.^2)+(1-m_u)./r_1+m_u./r_2);
contour(X, Y, P, [C_energy_target_orbit, C_energy_target_orbit], 'r')
xlabel('X')
ylabel('Y')
title('Earth parking orbit, Halo orbit, Transfer Orbit & Hill-Regions in 2-D')

plot(x(:,1), x(:,2), 'k')    
plot(x_halo(:,1), x_halo(:,2),'k','linewidth',1)
plot(h(3,1), h(3,2),'bx', h(1,1), h(1,2),'bx', h(2,1), h(2,2),'bx')
hold on 
plot(h(4,1), h(4,2),'bx', h(5,1), h(5,2),'bx')
plot(-m_u,0, 'go', 1-m_u, 0, 'go')

figure(2)
plot3(Xp(:,1)/384400000, Xp(:,2)/384400000, Xp(:,3)/384400000,'k','linewidth',1);
hold on
plot3(h(3,1), h(3,2) , 0, 'bx', h(1,1), h(1,2), 0, 'bx', h(2,1), h(2,2), 0, 'bx')
plot3(h(4,1), h(4,2), 0,'bx', h(5,1), h(5,2), 0,'bx')
plot3(-m_u,0, 0, 'go', 1-m_u, 0, 0,'go')
plot3(x(:,1), x(:,2), x(:,3), 'k')
plot3(x_halo(:,1), x_halo(:,2), x_halo(:,3), 'b')
axis([-1.2 1.2 -1.2 1.2 -0.006 0.006])
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Earth parking orbit, Halo orbit & Transfer Orbit in 3-D')

%Functions Used
%Jacobi energy
function C_energy = Jacobi_Const(y, v, mu)
r_1 = sqrt((mu+y(1,1))^2+(y(2,1))^2+(y(3,1))^2);
r_2 = sqrt((y(1,1)-(1-mu))^2+(y(2,1))^2+(y(3,1))^2);
C_energy = -(v(1,1)^2 + v(2,1)^2+v(3,1)^2)/2 + 2*((y(1,1)^2 + y(2,1)^2)/2 + (1-mu)/r_1 + mu/r_2);
end

function Xdot = satdyn(~,X)
mu_sun = 1.327e11;
r = X(1:3);
v = X(4:6);
Xdot(1:3) = v;
rn = norm(r);
Xdot(4:6) = -mu_sun/rn^3*r;
Xdot = Xdot(:);
return;
end

function y_dot = CRthreeBP(t, y, options, flag, G, mu)
r_1 = sqrt((mu+y(1))^2+(y(2))^2+(y(3))^2);
r_2 = sqrt((1-mu-y(1))^2+(y(2))^2+(y(3))^2);
m_1 = 1 - mu;
m_2 = mu;

y_dot=[y(4); 
       y(5); 
       y(6); 
       y(1)+2*y(5)+G*m_1*(-mu-y(1))/(r_1^3)+G*m_2*(1-mu-y(1))/(r_2^3); 
       y(2)-2*y(4)-G*m_1*(y(2))/(r_1^3)-G*m_2*y(2)/(r_2^3); 
       -G*m_1*y(3)/(r_1^3)-G*m_2*y(3)/(r_2^3)];
end

%STM at time t_final
function A = STM_2_CRTBP(t_0, t_f, x, mu)
t_span=[t_0,t_f];                       
y_0 = 0;
I_matrix = eye(6);         %Creating a 6*6 identity matrix                       
for n = 1:6
   for l = 1:6
      y(6*(n-1)+l) = I_matrix(n,l); 
   end    
end

y(37:42) = x';
y_0 = y';                                                    
options = odeset('RelTol',1e-13,'AbsTol',1e-22); 
[t,Y] = ode113('STM_CRTBP',t_span,y_0,options,[], mu);   %Integration

a = size(Y);        
b = a(1,1);         
c = Y(b,1:36);     %Finding the row number of t_final

d = 0;
for i = 1:6
    for j = 1:6
        d(i,j) = c(6*(i-1)+j);
    end    
end

A = d; %STM at t_final
end

function y_dot = STM_CRTBP(t, y, options, flag, mu)
G = 1;
x(1:3) = y(37:39);
Matrix = Matrix_G(x, mu);
Matrix_0 = zeros(3);
Matrix_I = eye(3);
K = [0, 2, 0; -2, 0, 0; 0, 0, 0];
derv_f = 0;      %Derivative of the N-body vector field f
derv_f = [Matrix_0, Matrix_I; Matrix, K];

for n = 1:6
    for l = 1:6
        A(n,l) = y(6*(n-1)+l);
    end
end

dervf_A = 0;
dervf_A = derv_f * A;

b = 0;
for i = 1:6
    for d = 1:6
        b(6*(i-1)+d) = dervf_A(i,d);
    end
end    

r_1 = sqrt((mu+y(37))^2+(y(38))^2+(y(39))^2);
r_2 = sqrt((1-mu-y(37))^2+(y(38))^2+(y(39))^2);
m_1 = 1 - mu;
m_2 = mu;

vector_field = [y(40); 
                y(41); 
                y(42);
                y(37)+2*y(41)+G*m_1*(-mu-y(37))/(r_1^3)+G*m_2*(1-mu-y(37))/(r_2^3); 
                y(38)-2*y(40)-G*m_1*(y(38))/(r_1^3)-G*m_2*y(38)/(r_2^3);
                -G*m_1*y(39)/(r_1^3)-G*m_2*y(39)/(r_2^3)];

y_dot = [b'; vector_field];
end

function f = CRTBP_vec_field(y, G, mu)
r_1 = sqrt((mu+y(1))^2+(y(2))^2+(y(3))^2);
r_2 = sqrt((y(1)-(1-mu))^2+(y(2))^2+(y(3))^2);
m_1 = 1 - mu;
m_2 = mu;

f=[y(4); 
    y(5); 
    y(6); 
    y(1)+2*y(5)+G*m_1*(-mu-y(1))/(r_1^3)+G*m_2*(1-mu-y(1))/(r_2^3); 
    y(2)-2*y(4)-G*m_1*(y(2))/(r_1^3)-G*m_2*y(2)/(r_2^3); 
    -G*m_1*y(3)/(r_1^3)-G*m_2*y(3)/(r_2^3)];
end

function h = Lagrange_Points_2(mu)
syms x_lp
mu = 0.012277471;
p = x_lp-((1-mu)*(x_lp+mu)/(abs(x_lp+mu)).^3)-(mu*(x_lp-1+mu)/(abs(x_lp-1+mu)).^3)==0;
x_lp = [vpasolve(p,x_lp,[0 1-mu]);vpasolve(p,x_lp,[1 1.2716]);vpasolve(p,x_lp,[-1.1984 -1]);0.5-mu;0.5-mu];
y_lp = [0;0;0;sqrt(3)/2;-sqrt(3)/2];

r_1 = zeros(1);
r_2 = zeros(1);
U = zeros(1);
C_J = zeros(1);

for i=1:5
    r_1(i) = sqrt((x_lp(i)+mu)^2+y_lp(i)^2);
    r_2(i) = sqrt((x_lp(i)-1+mu)^2+y_lp(i)^2);
    U(i) = -0.5*(x_lp(i)^2+y_lp(i)^2)-(1-mu)/r_1(i)-mu/r_2(i)-0.5*mu*(1-mu);
    C_J(i) = -2*U(i);
end
C_J;
h = [x_lp,y_lp];
end

