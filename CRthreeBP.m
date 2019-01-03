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