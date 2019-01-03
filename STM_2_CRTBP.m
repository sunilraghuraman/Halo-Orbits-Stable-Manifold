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
