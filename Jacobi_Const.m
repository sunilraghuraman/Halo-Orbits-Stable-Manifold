%Jacobi energy
function C_energy = Jacobi_Const(y, v, mu)
r_1 = sqrt((mu+y(1,1))^2+(y(2,1))^2+(y(3,1))^2);
r_2 = sqrt((y(1,1)-(1-mu))^2+(y(2,1))^2+(y(3,1))^2);
C_energy = -(v(1,1)^2 + v(2,1)^2+v(3,1)^2)/2 + 2*((y(1,1)^2 + y(2,1)^2)/2 + (1-mu)/r_1 + mu/r_2);

