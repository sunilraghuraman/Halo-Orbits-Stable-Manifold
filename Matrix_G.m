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
     
         
         
         
         
         