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