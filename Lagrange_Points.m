function Lagrange_Points(mu)
syms x
mu = 0.012277471;
p = x-((1-mu)*(x+mu)/(abs(x+mu)).^3)-(mu*(x-1+mu)/(abs(x-1+mu)).^3)==0;
x = [vpasolve(p,x,[0 1-mu]);vpasolve(p,x,[1 1.2716]);vpasolve(p,x,[-1.1984 -1]);0.5-mu;0.5-mu]
y = [0;0;0;sqrt(3)/2;-sqrt(3)/2]

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

C_J

plot(x(1),y(1),'r*',x(2),y(2),'r*',x(3),y(3),'r*',x(4),y(4),'r*',x(5),y(5),'r*')
hold on
xL = xlim;
yL = ylim;
line([0 0],yL);
line(xL,[0 0]);
title('Locations of all Lagrange points')
xlabel('x')
ylabel('y')

labels_1 = {'(0.83691512,0) L_1 ','(-1.00506264,0) L_3 ','(0.48784941,0.86602540) L_4 ','(0.48784941,-0.86602540) L_5 '};
x_plot_1 = [x(1),x(3),x(4),x(5)];
y_plot_1 = [0,0,sqrt(3)/2,-sqrt(3)/2];
text(x_plot_1, y_plot_1, labels_1,'VerticalAlignment','bottom','HorizontalAlignment','right');

text(x(2), 0, 'L_2 (1.15568216,0)','VerticalAlignment','bottom','HorizontalAlignment','left');
