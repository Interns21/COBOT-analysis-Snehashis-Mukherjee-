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

syms Pc1 Pc2 Pc3 Pc4 Pc5 real

%COM of each link of UR3e in m
Pc1 = [ 0; 0; 35.425]*10^(-3);
Pc2 = [ 64.70; 0; 0]*10^(-3);
Pc3 = [64.991; 0;0 ]*10^(-3);
Pc4 = [ 0; 9.368; -5.453]*10^(-3);
Pc5 = [0;-16.34;-1.8]*10^(-3);
Pc6 = [0;0;-1.159]*10^(-3);

syms m1 m2 m3 m4 m5 m6 I1 I2 I3 I4 I5 I6 th1ddot th2ddot th3ddot th4ddot th5ddot th6ddot g real

%angular acceleartion 
al0 = [0;0;0];
al1 = R01'*(al0 + cross(w0,[0;0;th1dot])) + [0;0;th1ddot];
al2 = R12'*(al1 + cross(w1,[0;0;th1dot])) + [0;0;th2ddot];
al3 = R23'*(al2 + cross(w2,[0;0;th1dot])) + [0;0;th3ddot];
al4 = R34'*(al3 + cross(w3,[0;0;th1dot])) + [0;0;th4ddot];
al5 = R45'*(al4 + cross(w4,[0;0;th1dot])) + [0;0;th5ddot];
al6 = R56'*(al5 + cross(w5,[0;0;th1dot])) + [0;0;th6ddot];

%linear acceleration
a0 = [0;g;0];
a1 = R01'*(a0 + cross(al0,P01) + cross(w0, P01));
a2 = R01'*(a1 + cross(al1,P12) + cross(w1, P12));
a3 = R01'*(a2 + cross(al2,P23) + cross(w2, P23));
a4 = R01'*(a3 + cross(al3,P34) + cross(w3, P34));
a5 = R01'*(a4 + cross(al4,P45) + cross(w4, P45));
a6 = R01'*(a5 + cross(al5,P56) + cross(w5, P56));

%acceleration of the COM
ac1 = a1 + cross(al1, Pc1) + cross(w1, cross(w1, Pc1));
ac2 = a2 + cross(al2, Pc2) + cross(w2, cross(w2, Pc2));
ac3 = a3 + cross(al3, Pc3) + cross(w3, cross(w3, Pc3));
ac4 = a4 + cross(al4, Pc4) + cross(w4, cross(w4, Pc4));
ac5 = a5 + cross(al5, Pc5) + cross(w5, cross(w5, Pc5));
ac6 = a6 + cross(al6, Pc6) + cross(w6, cross(w6, Pc6));

%moment of inertia of the links in kg-m^2
I1 = [1.22, 0, 0; 
      0, 1.259, 0;
      0, 0, 1.116]*10^(-3);
I2 = [.263, 0, 0;
      0, 1.512, 0;
      0, 0, 1.512]*10^(-3);
I3 = [.8698, 0, 0;
      0, 7.327, 0;
      0, 0, 7.327]*10^(-4);
I4 = [1.3, 0, 0;
      0, 1.039, .1653;
      0, .1653, 1.364]*10^(-4);
I5 = [9.82, 0, 0;
      0, 9.82, 0;
      0, 0, 1.271]*10^(-4);
I6 = eye(3)*10^(-4);

%Link Force and Moment in N and N-m
F1 = m1*ac1;
F2 = m2*ac2;
F3 = m3*ac3;
F4 = m4*ac4;
F5 = m5*ac5;
F6 = m6*ac6;

N1 = I1*al1 + cross(w1, I1*w1);
N2 = I2*al2 + cross(w2, I2*w2);
N3 = I3*al3 + cross(w3, I3*w3);
N4 = I1*al4 + cross(w4, I4*w4);
N5 = I1*al5 + cross(w5, I5*w5);
N6 = I1*al6 + cross(w6, I6*w6);

%Joint Force and Moment in N and N-m
f6 = [0;0;0];
f5 = R56*f6 + F5;
f4 = R56*f6 + F5;
f3 = R56*f6 + F5;
f2 = R56*f6 + F5;
f1 = R56*f6 + F5;

n6 = [0;0;0];
n5 = R56*n6 + cross(Pc5, F5) + cross(P56, (R56*f5)) + N5;
n4 = R45*n5 + cross(Pc4, F4) + cross(P45, (R45*f4)) + N4;
n3 = R34*n4 + cross(Pc3, F3) + cross(P34, (R34*f3)) + N3;
n2 = R23*n3 + cross(Pc2, F2) + cross(P23, (R23*f2)) + N2;
n1 = R12*n2 + cross(Pc1, F1) + cross(P12, (R12*f1)) + N1;

%tau expression of the joints
tau1 = simplify(n1(3));
tau2 = simplify(n2(3));
tau3 = simplify(n3(3));
tau4 = simplify(n4(3));
tau5 = simplify(n5(3));
tau6 = n6(3);

%inertia matrix
M = equationsToMatrix([tau1; tau2; tau3; tau4; tau5; tau6], [th1ddot, th2ddot, th3ddot, th4ddot, th5ddot, th6ddot]);
%gravitational matrix
G = equationsToMatrix([tau1; tau2; tau3; tau4; tau5; tau6], [g]);

