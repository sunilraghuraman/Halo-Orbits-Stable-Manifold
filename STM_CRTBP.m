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
