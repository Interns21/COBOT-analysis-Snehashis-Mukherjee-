close all
clc

syms alpha a theta d

syms th1 th2 th3 th4 th5 th6 real


N=6;

%Modified DH parameter table for UR3e
DHTABLE = [0,0,th1,.15185;
           pi/2,0,th2,0;
           0,-.2435,th3,0;
           0,-.2132,th4,.13105;
           -(pi/2),0,th5,.08535;
           (pi/2),0,th6,.0921;];

       
%Transfer Matrix
TDH = [cos(theta)       -sin(theta)*cos(alpha)       sin(theta)*sin(alpha)        a*cos(theta)
       sin(theta)       cos(theta)*cos(alpha)        -cos(theta)*sin(alpha)       a*sin(theta)
       0                sin(alpha)                   cos(alpha)                   d
       0                0                            0                            1];
   
A = cell(1:N);

for i=1:N
    alpha = DHTABLE(i,1);
    a     = DHTABLE(i,2);
    theta = DHTABLE(i,3);
    d     = DHTABLE(i,4);
    A{i}  = subs(TDH);
end

T = eye(4);

for i = 1:N
    T = T*A{i};
    T = simplify(T);
end

%Transfer matrixes for each coordinate system 
T01 = A{1};
T12 = A{2};
T23 = A{3};
T34 = A{4};
T45 = A{5};
T56 = A{6};

%position, normal, slide and approach vectors of end effector wrt base
%coordinate fame
p = simplify(T(1:3,4));
n = T(1:3,1);
s = T(1:3,2);
a = T(1:3,3);

%Rotational and Positional vectors 
R01 = A{1}(1:3,1:3);
P01 = A{1}(1:3,4);

R12 = A{2}(1:3,1:3);
P12 = A{2}(1:3,4);

R23 = A{3}(1:3,1:3);
P23 = A{3}(1:3,4);

R34 = A{4}(1:3,1:3);
P34 = A{4}(1:3,4);

R45 = A{5}(1:3,1:3);
P45 = A{5}(1:3,4);

R56 = A{6}(1:3,1:3);
P56 = A{6}(1:3,4);

syms th1dot th2dot th3dot th4dot th5dot th6dot real

%angular velocity
w0 = [0;0;0;];
w1 = R01'*(w0) + [0;0;th1dot];
w2 = R12'*(w1) + [0;0;th2dot];
w3 = R23'*(w2) + [0;0;th3dot];
w4 = R34'*(w3) + [0;0;th4dot];
w5 = R45'*(w4) + [0;0;th5dot];
w6 = R56'*(w5) + [0;0;th6dot];

w06 = simplify(R01*R12*R23*R34*R45*R56*w6);

%Linear Velocity
v0 = [0;0;0];
v1 = R01'*(v0+cross(w0,P01));
v2 = R12'*(v1+cross(w1,P12));
v3 = R23'*(v2+cross(w2,P23));
v4 = R34'*(v3+cross(w3,P34));
v5 = R45'*(v4+cross(w4,P45));
v6 = R56'*(v5+cross(w5,P56));

v06 = simplify(R01*R12*R23*R34*R45*R56*v6);

%Jacobian of the velocity equation
J = simplify(equationsToMatrix([v06,w06],[th1dot;th2dot;th3dot;th4dot;th5dot;th6dot]));