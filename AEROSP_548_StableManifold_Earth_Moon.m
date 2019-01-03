%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AEROSP 548 
%Final Project
%Team: Sunil Raghuraman, Monica A. Salunkhe
%Code Description: This code finds the stable manifolds when an initial 
%                  position of a periodic orbit is given. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
format long;

epsi = 2*10^(-8);
halo_period = 2.69239959528586;
T = 2.7*halo_period;
k = 20;
G = 1;
mu = 0.012277471;
r_0 = [0.83946302646687; 0; 0];
v_0 = [0; -0.02596831282986; 0];
Lagrange_Points_1;

t_0 = 0;
n_steps = 300;
t_span = linspace(0, halo_period+0.1*halo_period, n_steps);
y_0 = [r_0; v_0]';                                  
options = odeset('RelTol',1e-13,'AbsTol',1e-22);  
[t, x_halo] = ode45('CRthreeBP',t_span,y_0,options,[],G,mu);
Phi = STM_2_CRTBP(t_0, halo_period, [r_0; v_0], mu);  %fundamental matrix
[Phi_eigen_vectors, Phi_eigen_values] = eigs(Phi);

param_pts = round((n_steps-1)/k); %Getting parameterization 
                                  %points of the periodic orbit
for a = 1:k
halo_pt(a,1:6) = x_halo((a-1)*param_pts+1,1:6);
halo_t(a) = t((a-1)*param_pts+1);
end

I = eye(6);
unstable_eig = Phi_eigen_vectors(1:6,6);   %unstable eigenvector 
x_p_unstable(1,1:6) = halo_pt(1,1:6)'+epsi*I*unstable_eig; 

%Stable Orbit 
t_0 = 0;
n_steps = 2000;
t_span = linspace(0, T,n_steps);
y_0 = x_p_unstable(1,:)';                                     
options = odeset('RelTol',2.5e-13,'AbsTol',1e-22);    
[t,x_p_unstable1] = ode45('CRthreeBP_backward',t_span,y_0,options,[],G,mu);

for b = 2:k
    ini_conditions = b;
    perturbation_vec = STM_2_CRTBP(t_0, halo_t(b), [r_0; v_0], mu)*unstable_eig;
    perturbation_vec_cap = perturbation_vec/norm(perturbation_vec);
    x_p_unstable(b,1:6) = halo_pt(b,1:6)'+epsi*perturbation_vec_cap;
end

figure(1)
for j = 2:k
    manifoldLoop = j;
    t_0 = 0;
    n_steps = 2000;
    t_span = linspace(0, T, n_steps);
    y_0 = x_p_unstable(j,1:6)';                                     
    options = odeset('RelTol',2.5e-13,'AbsTol',1e-22);     
    [t,x_p_unstable_pos] = ode45('CRthreeBP_backward',t_span,y_0,options,[],G,mu);
    plot3(x_p_unstable_pos(:,1), x_p_unstable_pos(:,2), x_p_unstable_pos(:,3), 'b')
end
hold on 
plot3(1-mu, 0, 0,'ko')
plot3(-mu, 0, 0, 'ko')
plot3(x_halo(:,1), x_halo(:,2), x_halo(:,3), 'k.')


%Functions Used
function Lagrange_Points_1(mu)
syms x
mu = 0.012277471;
p = x-((1-mu)*(x+mu)/(abs(x+mu)).^3)-(mu*(x-1+mu)/(abs(x-1+mu)).^3)==0;
x = [vpasolve(p,x,[0 1-mu]);vpasolve(p,x,[1 1.2716]);vpasolve(p,x,[-1.1984 -1]);0.5-mu;0.5-mu];
y = [0;0;0;sqrt(3)/2;-sqrt(3)/2];

r_1 = zeros(1);
r_2 = zeros(1);
U = zeros(1);
C_J = zeros(1);

for i=1:5
    r_1(i) = sqrt((x(i)+mu)^2+y(i)^2);
    r_2(i) = sqrt((x(i)-1+mu)^2+y(i)^2);
    U(i) = -0.5*(x(i)^2+y(i)^2)-(1-mu)/r_1(i)-mu/r_2(i)-0.5*mu*(1-mu);
    C_J(i) = -2*U(i);
end

C_J;

plot(x(1),y(1),'r*',x(2),y(2),'r*',x(3),y(3),'r*',x(4),y(4),'r*',x(5),y(5),'r*')
hold on
xL = xlim;
yL = ylim;
line([0 0],yL);
line(xL,[0 0]);
title('Halo Orbit Stable Manifold Earth-Moon')
xlabel('x')
ylabel('y')

labels_1 = {'(0.83691512,0) L_1 ','(-1.00506264,0) L_3 ','(0.48784941,0.86602540) L_4 ','(0.48784941,-0.86602540) L_5 '};
x_plot_1 = [x(1),x(3),x(4),x(5)];
y_plot_1 = [0,0,sqrt(3)/2,-sqrt(3)/2];
text(x_plot_1, y_plot_1, labels_1,'VerticalAlignment','bottom','HorizontalAlignment','right');

text(x(2), 0, 'L_2 (1.15568216,0)','VerticalAlignment','bottom','HorizontalAlignment','left');
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

A = d;            %STM at t_final
end

function y_dot = STM_CRTBP(t, y, options, flag, mu)
G = 1;
x(1:3) = y(37:39);
Matrix = Matrix_G(x, mu);
Matrix_0 = zeros(3);
Matrix_I = eye(3);
K = [0, 2, 0; -2, 0, 0; 0, 0, 0];
derv_f = 0;       %Derivative of the N-body vector field f
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

function Matrix = Matrix_G(x, mu)
r_1 = sqrt((x(1)+mu)^2+x(2)^2+x(3)^2);
r_2 = sqrt((x(1)-(1-mu))^2+x(2)^2+x(3)^2);

%Derivative of the gradient of potential
U_1_x = 1-(1-mu)*(1/(r_1^3)-3*((x(1)+mu)^2)/(r_1^5))-mu*(1/(r_2^3)-3*((x(1)-(1-mu))^2)/(r_2^5));
U_1_y = 3*(1-mu)*x(2)*(x(1)+mu)/r_1^5+3*mu*x(2)*(x(1)-(1-mu))/r_2^5;
U_1_z = 3*(1-mu)*x(3)*(x(1)+mu)/r_1^5+3*mu*x(3)*(x(1)-(1-mu))/r_2^5;
U_2_x = U_1_y;
U_2_y = 1-(1-mu)*(1/(r_1)^3-3*x(2)^2/r_1^5)-mu*(1/r_2^3-3*x(2)^2/r_2^5);
U_2_z = 3*(1-mu)*x(2)*x(3)/r_1^5+3*mu*x(2)*x(3)/r_2^5;
U_3_x = U_1_z;
U_3_y = U_2_z;
U_3_z = (-1)*(1-mu)*(1/(r_1)^3-3*x(3)^2/r_1^5)-mu*(1/r_2^3-3*x(3)^2/r_2^5);
Matrix = [U_1_x, U_1_y, U_1_z; U_2_x, U_2_y, U_2_z; U_3_x, U_3_y, U_3_z];
end

function y_dot = CRthreeBP_backward(t, y, options, flag, G, mu)

r_1=sqrt((mu+y(1))^2+(y(2))^2+(y(3))^2);
r_2=sqrt((1-mu-y(1))^2+(y(2))^2+(y(3))^2);
m_1=1-mu;
m_2=mu;

y_dot=[-y(4); 
    -y(5); 
    -y(6); 
    -(y(1)+2*y(5)+G*m_1*(-mu-y(1))/(r_1^3)+G*m_2*(1-mu-y(1))/(r_2^3)); 
    -(y(2)-2*y(4)-G*m_1*(y(2))/(r_1^3)-G*m_2*y(2)/(r_2^3)); 
    -(-G*m_1*y(3)/(r_1^3)-G*m_2*y(3)/(r_2^3))];
end



