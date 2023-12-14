function [y]=mean_osculating_map(x,option,param)

% this code maps from mean to osculating element if option>0
% and osculating to mean if option<0 using Ref 2.
%
% [Ref. 2] Appendix G from, Junkins, John L., and Hanspeter Schaub. Analytical 
% mechanics of space systems. American Institute of Aeronautics and Astronautics, 2009.

a=x(1,1);
e=x(2,1);
inc=x(3,1);
bigO=x(4,1);
omega=x(5,1);
Mo=x(6,1);

mu = param.mu;
re = param.req;
J2 = param.j2;

if option>0
    % maps mean orbit elements 
    % to osculating orbit elements
    gam2=(J2/2)*(re/a)^2;
else
    %if the algorithm maps osculating orbit elements 
    % to mean orbit elements
    gam2=-(J2/2)*(re/a)^2;
end



if (Mo > -pi && Mo < 0) || Mo > pi
    E = Mo - e;
else
    E = Mo + e;
end

Epw = Mo;
for k =1:1000
    DeltaEpw = -(Mo-Epw+e*sin(Epw)) /(-1+e*cos(Epw));
    Epw = Epw + DeltaEpw;
    if DeltaEpw < 10e-14
        break
    end
end
E=Epw;

eta=sqrt(1-e^2);
gam2_p=gam2/eta^4;

theta=2*atan(sqrt((1+e)/(1-e))*tan(E/2));

ratio=(1+e*cos(theta))/(eta^2);


a_p=a+a*gam2*((3*cos(inc)^2-1)*((ratio)^3-(1/eta^3))+3*(1-cos(inc)^2)*(ratio)^3*cos(2*omega+2*theta));

Delta_e1=(gam2_p/8)*e*eta^2*(1-11*cos(inc)^2-40*((cos(inc)^4)/(1-5*cos(inc)^2)))*cos(2*omega);


term1=((3*cos(inc)^2-1)/(eta^6))*(e*eta+((e)/(1+eta))+3*cos(theta)+3*e*cos(theta)^2+e^2*cos(theta)^3);
term2=3*((1-cos(inc)^2)/(eta^6))*(e+3*cos(theta)+3*e*cos(theta)^2+e^2*cos(theta)^3)*(cos(2*omega+2*theta));
term3=gam2_p*(1-cos(theta)^2)*(3*cos(2*omega+theta)+cos(2*omega+3*theta));

Delta_e=Delta_e1+(eta^2/2)*(gam2*(term1+term2)-term3);

Delta_inc=-((e*Delta_e1)/(eta^2*tan(inc)))+(gam2_p/2)*cos(inc)*sqrt(1-cos(inc)^2)...
    *(3*cos(2*omega+2*theta)+3*e*cos(2*omega+theta)+e*cos(2*omega+3*theta));

angle_total=Mo+omega+bigO;
angle_total_p=angle_total+...
    (gam2_p/8)*eta^3*(1-11*cos(inc)^2-40*((cos(inc)^4)/(1-5*cos(inc)^2)))+...
    -(gam2_p/16)*(2+e^2-11*(2+3*e^2)*cos(inc)^2-40*(2+5*e^2)*((cos(inc)^4)/(1-5*cos(inc)^2))-400*e^2*((cos(inc)^6)/(1-5*cos(inc)^2)^2) )+...
    (gam2_p/4)*(-6*(1-5*cos(inc)^2)*(theta-Mo-e*sin(theta))+(3-5*cos(inc)^2)*(3*sin(2*omega+2*theta)+3*e*sin(2*omega+theta)+e*sin(2*omega+3*theta)) )+...
    -(gam2_p/8)*e^2*cos(inc)*(11+80*((cos(inc)^2)/(1-5*cos(inc)^2)) +200*((cos(inc)^4)/(1-5*cos(inc)^2)^2)  )+...
    -(gam2_p/2)*cos(inc)*(6*(theta-Mo+e*sin(theta))-3*sin(2*omega+2*theta)-3*e*sin(2*omega+theta)-e*sin(2*omega+3*theta) );

eDelta_M=(gam2_p/8)*e*eta^3*(1-11*cos(inc)^2-40*((cos(inc)^4)/(1-5*cos(inc)^2)) )+...
    -(gam2_p/4)*eta^3*(  2*(3*cos(inc)^2-1)*(ratio^2*eta^2+ratio+1)*sin(theta)+ ...
    3*(1-cos(inc)^2)*( (-ratio^2*eta^2-ratio+1)*sin(2*omega+theta)+...
    (ratio^2*eta^2+ratio+1/3)*sin(2*omega+3*theta)  ));


Delta_bigO=-(gam2_p/8)*e^2*cos(inc)*(11+80*(cos(inc)^2)/(1-5*cos(inc)^2)+200*(cos(inc)^4)/(1-5*cos(inc)^2)^2)+...
           -(gam2_p/2)*cos(inc)*(6*(theta-Mo+e*sin(theta)) -3*sin(2*omega+2*theta)-3*e*sin(2*omega+theta)-e*sin(2*omega+3*theta));


d1=(e+Delta_e)*sin(Mo)+(eDelta_M)*cos(Mo);
d2=(e+Delta_e)*cos(Mo)-(eDelta_M)*sin(Mo);

Mo_p=atan(d1/d2);

e_p=sqrt(d1^2+d2^2);

d3=(sin(inc/2)+cos(inc/2)*(Delta_inc/2))*sin(bigO)+sin(inc/2)*Delta_bigO*cos(bigO);
d4=(sin(inc/2)+cos(inc/2)*(Delta_inc/2))*cos(bigO)-sin(inc/2)*Delta_bigO*sin(bigO);

bigO_p=atan(d3/d4);

inc_p=2*asin(sqrt(d3^2+d4^2));

omega_p=angle_total_p-Mo_p-bigO_p;

y=[a_p;e_p;inc_p;bigO_p;omega_p;Mo_p];
