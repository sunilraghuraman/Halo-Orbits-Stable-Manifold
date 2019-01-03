clear 
format long



%--------------------------------------------------------------------------
%----------------------NOTES-----------------------------------------------
%--------------------------------------------------------------------------

%Dependencies:
%This program calls 'nBodyWpar.m', and 'sitnikowVewton.m'

%USE:
%Program simulates the Sitnikov Problem.
%This program has many userdefined parameters.  Several define the shape of
%the orbits of the primaries, which are ellipses.  Several more have to do
%with computing the poincare map of the motion of the third body.  They are
%marked clearly below but are located in TWO seprate sections.

%Once the user sets these (or leaves them as they are now)
%then save and type 'sitnikovMap' at the MatLab command prompt with the
%present program and it's dependencies in the search path.


%OUTPUT:
%The output is several plots.

%One plot shows the orbits of the three bodies in physical space.
%One shows the height of the third body abov or below the xy plane
%as a function of time.

%Then there are two Poincare maps.  One in angular coordinates.
%and the other lifted to the square.

%PURPOSE:
%The purpose of the program is to study the Poincare mapping associated
%with the Sitnikov problem.  The program was used to generate the results
%about this problem at the end of the first set on notes where the Sitnikov
%proplem is discussed in some detail.  The program lets you decide how many
%initial conditions to simulate, and how many intersections with the Poincare
%section to compute for each initial condition.



%------------------------DATA----------------------------------------------

N=3;        %number of bodies
%gravatational constant
G=1;
%Masses
m1=5.0;
m2=5.0;
m3=0.00000001;
Mass=[m1 m2 m3];

%--------------------------------------------------------------------------
%----------------------SHAPE PARAMETERS------------------------------------
%---------------------------FOR THE----------------------------------------
%-----------------------PRIMARY BODIES--------------------------------------
%--------------------------------------------------------------------------



%--------Compute the appropriate initial conditions for the primaries------

e=0.65           %eccentricity: Should be between zero and one (but not 1).
                 %              For e=0 the problem is completely
                 %              integrable.  e can be considered the
                 %              perturbation parameter for the problem.
d=5              %max distance
M=m1+m2

%---------------Compute the needed two body initial conditions-------------

%Kepler constants
Ec2=(G*M)^2*(e^2-1)/2
c2=d*G*M*(1-(1+(2*Ec2)/(G*M)^2)^(1/2))
c=sqrt(c2)
thetaDotZero=c/d^2
v_0=d*thetaDotZero

p=c2/(G*M)
a=p/(1-e^2)

%conversion to two body coordinates
x2_0=d/(1+m2/m1)
x1_0=-(m2/m1)*x2_0

x2Dot_0=v_0/(1+m2/m1)
x1Dot_0=-(m2/m1)*x2Dot_0


T=2*pi*a^(3/2)/sqrt(G*M)


%--------------------------------------------------------------------------




%----------------------CONTROL PARAMETERS----------------------------------
    %NOTE:  The first 4 variables are set by the user

    %This variable decides how many initial conditions will be used to
    %generate the Poincare map, or how many points in the domain of the
    %mapping to consider.
        k=300

    %This variable decides how many crossings of the xy-plane to keep for
    %every initial condition.  This is equivalent to specifying how many
    %itterates of the poincare map to keep for every point.
        iterates=1000

    %NOTE: To get really good looking graphics one should use k=100-300 and
    %poincareItterates=1000.  However the program will run for a very long
    %time with these values (days even).  As they are set now the program
    %runs in a few minutes and finds an invariant torus or two.

    %NOTE:  zPrime=0 is a fixed point of the map.  Setting the variable
    %below is deciding how far from the fixed point you want to map.
    %The boundary between orbits which escape to infinity and those
    %which are bounded is near z_0=3.643815.  After this orbits escape to
    %infinity.

    %initial z velocity
    zPrime_0=0.5

    %This determines the length of interval of initial velocities which
    %will be used
    zPrime_max=1.5



    %NOTE: The following parameters should not need to be changed.

    %sets the tolerence for the map.  The poincare map looks for crossings
    %of the xy plane.  epsilon determines how close a trajectory point must
    %be to the plane to be considered "on" the plan.  Standard value is
    %ten to the minus six.
        epsilon=0.000001

   
    
    %The program integrates for a time, then checks for crossings.  It repeats 
    %process untill it has found the desired number of crossings.  The variable
    %below defines how long to integrate before looking for crossings.  Thsi is
    %somewhat arbitrary.

        timeStep=1.0;

    %The poincare map is no a bounded map on all of R^2.  There is 
    %z_0' above which trajectories escape to infinity.  Then the domain of
    %the map is a topological disk whos boundary we do not know a priori.  
    %This has to be garded against.  The while loop will stop executing if 
    %the z component of some trajectory gets too big.  If this happens we will
    %conclude that we are out of domain and the program stops. 
             outOfDomain=0;
             escape=8;      %We consider the trajectory to have escaped if 
                             %the z or z' component exceeds this tolerence.
             boundDomain=0;  %We will also keep track of the last trajectory
                             %that was in the domain.  This gives an
                             %estimate as to the location of the boundary    

%--------------------------------------------------------------------------


%---------------Derived Parameters-----------------------------------------
%-----------(Should not need to touch)-------------------------------------



%period of the primaries
    P=T;

      
    
 %NOTE: The following parameters are necessary but depend on the above

    %if k=1 no steps are necessare and the following variable makes
    %no sense.
    if k>1
        %determine the step size between successive initial conditons
        zStep=(zPrime_max-zPrime_0)/(k-1);
    else
        %if k=1 just set step to zero so it is at least defined later
        zStep=0;
    end %end seeing if steps are necessary

%INITIALIZE:
    %This is the number of iterates of the current initial conditions
    %found this far through the while loop
    z_nIterates=0;
    %This is the total number of points in the map so far
    totalIterates=0;
    %Local name of initial conditions:
    initial_n=0;
    %keep track of how long an initial condition is integrated
    trajectoryTime=0;
    
    %number of times through the while loop
    checkWhile=0;
    %number o times the various loops are called
    interpTimes=0;
    ReIntegrate=0;
    checkTotal=0;
    calledNewton=0;

%--------------------------------------------------------------------------
%----------------------COMPUTATION OF POINCARE MAP-------------------------
%--------------------------------------------------------------------------

for n=1:k
    %I like to see what loop I'm in at command line, just to see how fast
    %the program is running and where it is.  Comment this out or delete it
    %if you prefer.  Same goes for any 'check' variable
    check_n=n
    checkTotal=totalIterates
    %First the position
    initial_n = [x1_0; 0.0; 0.0;
          x2_0; 0.0; 0.0;
          0.0;  0.0; 0.000000001;
          0.0 ; x1Dot_0; 0.0;
          0.0; x2Dot_0; 0.0;
          0.0; 0.0; zPrime_0+(n-1)*zStep]';
            
    %THE MAIN WHILE LOOP:  This computes the number of times the
    %trajectory with initial condition 'initial_n' crosses the poincare
    %section with positive y velocity. The desired number of such crossings
    %is given by the integer 'iterates'.  The count of how many have been
    %found so far is in the raviable:
              x_nIterates=0;   %This should be set back to zero before
                               %every run through the while loop.
   
    %set trajectory time back to zero for each new initial condition
             trajectoryTime=0;
    
   
               
    %Begining of while loop
    while (x_nIterates <=iterates) & (outOfDomain ~= 1)
        checkWhile=checkWhile+1;
               
        %Intigrate this for time 'timeLength'.
           tspan=[0 timeStep];
           options=odeset('RelTol',1e-12,'AbsTol',1e-12);
           [t,trajectory_n] = ode113('nBodyWpar',tspan,initial_n,options,flag,N,G,Mass);
           temp=size(trajectory_n);
           numSteps=temp(1,1);
           
                    
       if ((abs(trajectory_n(numSteps,9))>=2*escape) | (abs(trajectory_n(numSteps,18))>=escape))
           outOfDomain=1;
           boundDomain=trajectory_n(1,18);
       end %end of outOfdomain conditional
           
           
       %Now check for pairs of points along the trajectory between which
       %the sign of y changes and the sign of ydot is positive
       for i=2:numSteps
           
           if (outOfDomain == 1)
               break
           end    
           
           
           %Check for sign change of y and correct sign of ydot.
           if (sign(trajectory_n(i,9)) ~= sign(trajectory_n(i-1,9)))  
                %note the approximate time of the crossing
                roughTime=t(i-1);
                %It could happen that one of these two points is already a
                %crossing to the desired tolerence.  We check for this:
                if (abs(trajectory_n(i-1,9))<epsilon)
                    intersection=trajectory_n(i-1,:);
                    intersectionTime=t(i-1);
                elseif (abs(trajectory_n(i,9))<epsilon)
                    intersection=trajectory_n(i,:);
                    intersectionTime=t(i);
                else
                %If (as is most likely the case) neither of the points 
                %for which the sign change occurs are in tolerence
                %(to be considered as crossings), 
                %then a Newton method determines the 
                %crossing to the desired accuracy.
                    [timeCorrection, intersection]=sitnikovNewton(t(i-1),trajectory_n(i-1,:),t(i),trajectory_n(i,:),epsilon,N,G,Mass);
                    calledNewton=calledNewton+1;  %(Print to scren)
                end %end of the intersection finding conditional.
                   
                %Keep track of how many iterates of the current initial
                %condition have been found
                    x_nIterates=x_nIterates+1; 
                %Keep track of how many total intersections have been found
                    totalIterates=totalIterates+1 %(print to screne)
                %Stor the intersection information for plotting later
                    %The actual time of the crossing is
                    moserTime=trajectoryTime+roughTime+timeCorrection;
                    %store the poincarePoint
                    moser(totalIterates,1)=mod(2*pi*moserTime/P,2*pi);
                    moser(totalIterates,2)=abs(intersection(18));
           end %end of pair checking if statement
       end %end of pair checking loop (i=2:numSteps index)
                      
        %make sure variables are ready for next while loop: unless all the
        %iterates have been found we will go fhrought the while loop again.
        %In that case the new initial condition is the final condition of
        %trajectory_n
              
        %How long has this initial condition been integrated
        
        initial_n=trajectory_n(numSteps,:);
        trajectoryTime=trajectoryTime+timeStep;
    
    end %end of main while loop (x_nIterates <= iterates)
           
    %The set up for the next time through the 'n' loop is done at the
    %begining of the loop instead of here at the end.
    %So just end and go back.
          
        %if the trajectories are blowing up then break
        if (outOfDomain == 1)
            break
        end %end break conditional   
            
end  %end of outmost for loop (n=1:k index)


%--------------------------------------------------------------------------
%--------------------ANALYSIS OF PERFORMANCE-------------------------------
%--------------------------------------------------------------------------

%We want to know how many meaning full intersections have been computed.
%This can be found by subtracting from the number of points in
%'PoincareMap' the number of times the point [-mu, 0, 0, 0, 0, 0] occurs.
%This is described in the comments of the program 'CRTBPpoincareNewton'.
%Then count these now:

    %'PoincareMap' will contain 'k'*'iterates' intersections.
    computedIntersections = k*iterates
    %check this:
    Temp=size(moser)
    checkPMsize=Temp(1,1)
    checkTotal = totalIterates
    
    %All three of the above should be the same number. If so then the
    %variable sizesGood will be one.
    if ((computedIntersections+1)==checkPMsize) & (checkPMsize==checkTotal)
        sizeGood=1
    else 
        sizeGood=0
    end    
    
    %initialize the counter
    falsePositives=0

for p=1:computedIntersections
    if moser(p,:)==[0 0]
        falsePositives=falsePositives+1;
    end %end of falsePositives conditional loop    
end %end of computedIntersections for loop   

%output the results
checkFalsePositives=falsePositives
trueIntersections = computedIntersections - falsePositives


%How many times was the interpolator called?
checkInterpTimes=interpTimes

%How many times is the ReIntegrate loop being called?
checkReIntegrate=ReIntegrate

%bound domain: If this is zero then the boundary was not found
%If it is nonzero then it is the z' coordinate of the last orbit that did
%not escape.
checkBound=boundDomain


save sitnikovData1 moser  
    
    
    
    
    
    
    

%--------------------------------------------------------------------------

%If running on UNIX in Background the following should be commented out


%--------------------------------------------------------------------------
%-------------------------PLOTTING-----------------------------------------
%--------------------------------------------------------------------------


%--------Integrate the primarise and have a look at their orbits----------

%integrate
%numSteps = 4000;
%tf=moserTime;
%tspan = linspace(0, tf, numSteps);
%initial positions and velocities
%options=odeset('RelTol',1e-12,'AbsTol',1e-12);
%[t, y] = ode113('nBodyWpar',tspan,y0,options,flag,N,G,Mass);


%hold on
%plot3(y(:,1), y(:,2), y(:,3),'b', y(:,4), y(:,5), y(:,6), 'g',  y(:,7), y(:,8), y(:,9), 'r')
%title('physical configuration of the Sitnikov problem')

%figure
%hold on
%plot(t(:),y(:,9),'b')
%title('height of the third body vs. time')

%--------------------------------------------------------------------------



%------------------------Moser Mapping-------------------------------------
figure
plot(moser(:,1),moser(:,2),'b.')
title('lift of the Poincare map to the square')

figure
polar(moser(:,1),moser(:,2),'b.')
title('Poincare section of the Sitnikov problem')
%------------------------PROGRAM OVER--------------------------------------
