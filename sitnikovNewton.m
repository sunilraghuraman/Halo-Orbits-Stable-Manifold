function varargout=sitnikovNewton(t0, initialPoint, t1, finalPoint, tolerence, N,G,Mass)


%--------------------------------------------------------------------------
%-------------------------NOTES--------------------------------------------
%--------------------------------------------------------------------------

%INPUT:
%This function is passed data about two points, one on each side of the
%line y=0.  Both are points on a trajectory, and 'initialPoint' is the one
%that occurs first.  The time that each point occurs on the trajectory is
%passed into the functio, but only the difference will matter.  Finally, an
%error tolerence is passed in as well.  This tells the function how close a
%result has to be to the x axis so that the Newton method can be considered
%to have converged.  Ofcourse the mu and G being used should be passed.

%OUTPUT:
%The purpose of the function is to find a point on the trajectory between
%'initialPoint'and 'finalPoint' which is an intersection with the x axis
%up to the allowed tolerence.  This is the 'x' that the function returns.

%HOW IT WORKS:
%As the name suggeste the function implements a Newton method to find the
%point.  However since we want a new trajectory that begins at say
%'initialPoint' and flows to the x axis, we only want to find the time that
%'initialPoint' must be flowed.  The initila guess is t1-t0.  Then the
%Newton method is only one dimensional.  Since it will only need the
%derivative of the flow with respect to time (and not w.r.t. initial
%conditions) all that is required is an evaluation of the vector field (and
%there is no need to integrate the variational equations).


%A STICKY POINT:
%Of the two points 'initialPoint' and 'finalPoint' one will be closer to
%the x axis than the other (possibly much closer).  Then the progrma will
%pick the point that is closer, and use that as the initial guess for the
%crossing.  The subtly is that if 'initialPoint' is closer,  then we want
%to examine trajectories begining from 'finalPoint' and flowing back to
%'initialPoint', as the one we want is one that flows not quite to
%'initialPoint' but intersects the x axis.  In this case the equations of
%motion for the CRTBP are not integrated, but instead the equations for the
%reverse time problem.  Then the function must code two different Newton
%methods and decide which one to use.




%--------------------------------------------------------------------------
%----------------------Computation-----------------------------------------
%--------------------------------------------------------------------------

%Check while will count the number of times the Newton step excutes within
%the algorithm.  If the step excutes more than say 15 times, then the
%method is probably not going to converge.  In that case we do not want to
%be stuck in the Newton loop forever.  If 'checkWhile' gets bigger than 
%'breakWhile' then the Newton algorithm should break and return some 
%default value 'newtonFailed'.  If 'newtonFailed' is some value that the 
%Poincare map will never otherwise take on, then it gives a good count of 
%how many times the algorithm failed to converge.
checkWhile=0;
breakWhile=16;
maxTime=8;                            %If the Newton algorithm is not 
                                       %converging the integration times 
                                       %may get longer and longer.  However
                                       %under reasonable circumstances the
                                       %initial time guess should be small
                                       %and get smaller.  Then if t_n is
                                       %getting bigger than say 10, the
                                       %algorithm is almost deffinitiel
                                       %diverging.  Check for this as well.

newtonFailed=zeros(1, 18);     %This is the location of the primary 
                                       %which is a singularity of the
                                       %system.  There is no way an
                                       %itterate of point could ever reach
                                       %this value, so if we count the
                                       %number if times it is returned then
                                       %we get a good count of how many
                                       %times the algorithm failed to
                                       %converge.


%Regardless of which point we begin at the time it took to get from one to
%the other is the same.  This is the initial guess for the time from start
%to intersection.

deltaT=t1-t0;
t_n=t1-t0;

%Determine which point is closest and find the intersection:
if abs(finalPoint(9))<=abs(initialPoint(9))
    %In this case 'finalPoint' is the initial guess as to the location of
    %the crossing.  Then we begin with trajectories starting at
    %'initialPoint' and use the forward time dynamics with a Newton
    %algorithm to refine the guess.

    %The initial guess for the location of the crossing is:
    y_n=finalPoint(9);
      
    %The Forward Time Newton Algorithm: excitue untill tolerence is met
    while abs(y_n) >= tolerence
        if (checkWhile < breakWhile) & (t_n < maxTime)
        %Integrate.  fx_n=f_n(endTime,2) is f(x_n) in the Newton Algorithm
             tspan=[0 t_n];
             options=odeset('RelTol',1e-12,'AbsTol',1e-12);
             [t,f_n] = ode113('nBodyWpar',tspan,initialPoint,options,flag,N,G,Mass);
        %Compute 'endTinme':
             sizef_n=size(f_n);
             endTime=sizef_n(1,1);
        %Compute fx_n
             fx_n=f_n(endTime,9);
        %Compute Dfx_n = [f(x_n)']^-1 =f_n(endTime, 5):
             Dfx_n = f_n(endTime, 18);
        %Compute the Newton Step
             t_n = t_n - fx_n/Dfx_n;
        %This is the new time estimate.  Now compute the new crossing
        %estimate:
             %First Integrate over new time;
                  tspan=[0 t_n];
                  options=odeset('RelTol',1e-12,'AbsTol',1e-12);
                  [t,f_n] = ode113('nBodyWpar',tspan,initialPoint,options,flag,N,G,Mass);
            %Now y_n is the second component of the end of this:
                  %Recompute 'endTinme':
                  sizef_n=size(f_n);
                  endTime=sizef_n(1,1);
                  y_n=f_n(endTime,9);
                  crossing_n=f_n(endTime,:);
                  timeCorrection=t_n;
        %increment 'checkWhile'
        checkWhile=checkWhile+1;
        else
           crossing_n=newtonFailed;  
           y_n = 0.1*tolerence;       
           timeCorrection=0;
        end
    end %end of while: the Forward time Newton Algorithm
else
    %In this case 'initialPoint' is closer to the x-axis so it is our first
    %guess.  Then we flow from 'finalPoint' under the backwards time
    %dynamics, and use a Newton algorithm to refine the guess.
         
    %The initial guess for the location of the crossing is:
    y_n=initialPoint(9);

    %The Backward Time Newton Method: excute untill tollerence is met
    while abs(y_n) >= tolerence
        if (checkWhile < breakWhile) & (t_n < maxTime)
        %Integrate.  fx_n=f_n(endTime,2) is f(x_n) in the Newton Algorithm
             tspan=[0 t_n];
             options=odeset('RelTol',1e-12,'AbsTol',1e-12);
             [t,f_n] = ode113('nBodyWpar',tspan,finalPoint,options,flag,N,G,Mass);
        %Compute 'endTinme':
             sizef_n=size(f_n);
             endTime=sizef_n(1,1);
        %Compute fx_n
             fx_n=f_n(endTime,9);
        %Compute Dfx_n = [f(x_n)']^-1 =f_n(endTime, 5):
             Dfx_n = -f_n(endTime, 18);
        %Compute the Newton Step
             t_n = t_n - fx_n/Dfx_n;
        %This is the new time estimate.  Now compute the new crossing
        %estimate:
             %First Integrate over new time;
                 tspan=[0 t_n];
                 options=odeset('RelTol',1e-12,'AbsTol',1e-12);
                 [t,f_n] = ode113('nBodyWpar',tspan,finalPoint,options,flag,N,G,Mass);
            %Now y_n is the second component of the end of this:
                  %Reompute 'endTinme':
                  sizef_n=size(f_n);
                  endTime=sizef_n(1,1);
                  y_n=f_n(endTime,9);
                  crossing_n=f_n(endTime,:);
                  timeCorrection=deltaT-t_n;
            %increment 'checkWhile'
            checkWhile=checkWhile+1;
        else
             crossing_n=newtonFailed;  
             y_n=0.1*tolerence;
             timeCorrection=0;
        end  %end of check while if conditional    
    end %end of while: the backward time Newton Algorithm
end %end of the if then conditional

%Return the point of intersection to the caller


varargout(1,1)={timeCorrection};
varargout(1,2)={crossing_n};









