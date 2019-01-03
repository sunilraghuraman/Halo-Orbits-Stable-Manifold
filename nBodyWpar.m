function ydot=nBodyWpar(t,y,options,flag,N,G,Mass)



%The following program computes the right hand side of the N-Body Problem.
%In addition to the usual suspects you have to pass the numbe of bodies 'N'
%the value of the gravitational constant 'G', and a row vector containing
% the masses.  This vector is locally named 'Mass' (often this is a good
%global name for it as well).

%Enter the acceleration components (bottom half of RHS)
acc=0;                      %Initialize the variables
s1=0;
s2=0;
s3=0;
for i=1:N
    for j=1:N
            
Rij=(y(3*i-2)-y(3*j-2))^2+(y(3*i-1)-y(3*j-1))^2+(y(3*i)-y(3*j))^2;
            %compute the three components of acceleration i
            if j~=i
               s1=s1+(Mass(1,j)/(sqrt(Rij))^3)*(y(3*j-2)-y(3*i-2));
               s2=s2+(Mass(1,j)/(sqrt(Rij))^3)*(y(3*j-1)-y(3*i-1));
               s3=s3+(Mass(1,j)/(sqrt(Rij))^3)*(y(3*j)-y(3*i));
            else
               s1=s1+0;
               s2=s2+0;
               s3=s3+0;
            end
    end
    acc(3*i-2)=G*s1;
    acc(3*i-1)=G*s2;
    acc(3*i)=G*s3;
    s1=0;
    s2=0;
    s3=0;
end

%Store the accelerations
accelerations=0;
accelerations=acc;

%enter the velociety components
velocities=0;
for i=1:(3*N)
velocities(i)=y(3*N+i);
end

%constructs a vector whose first 3*N entries are the velocities and whose
%last 3*N entries are the accelerations
Yprime=0;
Yprime(1:3*N,1)=velocities;
Yprime(3*N+1:6*N,1)=accelerations;

%Hands the right hand side to the integrator
ydot=Yprime;

