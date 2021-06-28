close all
clc

syms nx ny nz sx sy sz ax ay az Px Py Pz

syms th1 th2 th3 th4 th5 th6 real

%DH parameter table
DHTABLE = [pi/2,0,th1,.15185;
           0,-.2435,th2,0;
           0,-.2132,th3,0;
           pi/2,0,th4,.13105;
           -(pi/2),0,th5,.08535;
           0,0,th6,.0921;];

%transfer matrix    
TDH = [0,0,-1,779/20000;
       0,1,0,0;
       1,0,0,181/5000;
       0,0,0,1;];

%position of the 5th joint 
Px5 = TDH(1,4) - DHTABLE(6,4)*TDH(1,3);
Py5 = TDH(2,4) - DHTABLE(6,4)*TDH(2,3);

s1 = [sqrt(Px5^2 + Py5^2 - DHTABLE(4,4)^2), -sqrt(Px5^2 + Py5^2 - DHTABLE(4,4)^2)];

%calculation of theta1 value
theta1 = atan2(Px5, -Py5) + atan2(s1(1),DHTABLE(4,4));                     

m = TDH(1,4)*sin(theta1) - TDH(2,4)*cos(theta1);
m1 = (m - DHTABLE(4,4))/DHTABLE(6,4);

%calculation of theta5 value
theta5 = m1*atan2(sqrt(1-m1^2),m1);

m2 = TDH(3,3)/sin(theta5);

theta234 = atan2(m2,sqrt(1-m^2));

if isnan(theta234)
   disp("The angles theta2 theta3 theta4 theta6 are  undefined and the joints are parallel to each other")
else
    m3 = -sin(theta234)*cos(theta1);
    s2 = [sqrt(TDH(1,1)^2+TDH(1,2)-m3^2),-sqrt(TDH(1,1)^2+TDH(1,2)-m3^2)];
    
    %calculation of theta6 value
    theta6 = atan2(TDH(1,1),TDH(1,2)) + atan2(s2,m3);
    m4 = -(TDH(1,4)*cos(theta1) + TDH(2,4)*sin(theta1) + DHTABLE(6,4)*cos(theta234)*sin(theta5) - DHTABLE(5,4)*sin(theta234));
    m5 = -(TDH(3,4) - DHTABLE(1,4) + DHTABLE(6,4)*sin(theta5) + DHTABLE(4,4)*cos(theta234));
    s3 = (m4^2 + m5^2 - .2132^2 - .2435^2)/(2*.2132*.2435);
    
    %calculation of theta3 value
    theta3 = atan2(sqrt(1-s3^2),s3);
    
    %calculation of theta2 value
    theta2 = atan2(((DHTABLE(2,2)*cos(theta3)+DHTABLE(3,2))*m5 - DHTABLE(2,2)*sin(theta3)*m4),((DHTABLE(2,2)*cos(theta3)+DHTABLE(3,2))*m4 - DHTABLE(2,2)*sin(theta3)*m5));
    
    %calculation of theta4 value
    theta4 = theta234 -theta2 - theta3;
end