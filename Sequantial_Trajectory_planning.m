close all
clc

%initial and final position conditions
xo = .050; yo = .038; zo = .055;
xf = .055; yf = .043; zf = .050;

%initial and final velocity condition
xdoto = 0; ydoto = 0; zdoto = 0; wxo = 0; wyo = 0; wzo = 0; 
xdotf = 0; ydotf = 0; zdotf = 0; wxf = 0; wyf = 0; wzf = 0; 

%time duration
tt = 18;
tf = 18/6;
t = 1:.1:tt;

%theta initial and final value
thetao = ThetaValues(0,0,1,0,1,0,-1,0,0,xo, yo, zo);
thetaf = ThetaValues(0,0,1,0,1,0,-1,0,0,xf, yf, zf);

%intial and final jacobian
Jo = Jacobian(thetao(1), thetao(2), thetao(3), thetao(4), thetao(5), thetao(6));
Jf = Jacobian(thetaf(1), thetaf(2), thetaf(3), thetaf(4), thetaf(5), thetaf(6));

%initial and final thetadot values
qdoto = Jo\[xdoto; ydoto; zdoto; wxo; wyo; wzo];
qdotf = Jf\[xdotf; ydotf; zdotf; wxf; wyf; wzf];

%cubic ploynomial matrix
A = [1, 0, 0, 0;
     0, 1, 0, 0;
     1, tf, tf^2, tf^3;
     0, 1, 2*tf, 3*tf^2];
 
b1 = [thetao(1); qdoto(1); thetaf(1); qdotf(1)];
b2 = [thetao(2); qdoto(2); thetaf(2); qdotf(2)];
b3 = [thetao(3); qdoto(3); thetaf(3); qdotf(3)];
b4 = [thetao(4); qdoto(4); thetaf(4); qdotf(4)];
b5 = [thetao(5); qdoto(5); thetaf(5); qdotf(5)];
b6 = [thetao(6); qdoto(6); thetaf(6); qdotf(6)];

tc1 = A\b1;
tc2 = A\b2;
tc3 = A\b3;
tc4 = A\b4;
tc5 = A\b5;
tc6 = A\b6;

for i = 1:length(t)
    if t(i) <= tf
        th1(i) = [1, t(i), t(i)^2, t(i)^3]*tc1;
        th2(i) = thetao(2);
        th3(i) = thetao(3);
        th4(i) = thetao(4);
        th5(i) = thetao(5);
        th6(i) = thetao(6);
        
    elseif t(i) <= 2*tf
        th1(i) = thetaf(1);
        th2(i) = [1, t(i)-tf, (t(i)-tf)^2, (t(i)-tf)^3]*tc2;
        th3(i) = thetao(3);
        th4(i) = thetao(4);
        th5(i) = thetao(5);
        th6(i) = thetao(6);
    
    elseif t(i) <= 3*tf
        th1(i) = thetaf(1);
        th2(i) = thetaf(2);
        th3(i) = [1, t(i)-2*tf, (t(i)-2*tf)^2, (t(i)-2*tf)^3]*tc3;
        th4(i) = thetao(4);
        th5(i) = thetao(5);
        th6(i) = thetao(6);
    
    elseif t(i) <= 4*tf
        th1(i) = thetaf(1);
        th2(i) = thetaf(2);
        th3(i) = thetaf(3);
        th4(i) = [1, t(i)-3*tf, (t(i)-3*tf)^2, (t(i)-3*tf)^3]*tc4;
        th5(i) = thetao(5);
        th6(i) = thetao(6);
        
      elseif t(i) <= 5*tf
        th1(i) = thetaf(1);
        th2(i) = thetaf(2);
        th3(i) = thetaf(3);
        th4(i) = thetaf(4);
        th5(i) = [1, t(i)-4*tf, (t(i)-4*tf)^2, (t(i)-4*tf)^3]*tc5;
        th6(i) = thetao(6);
        
       elseif t(i) <= 6*tf
        th1(i) = thetaf(1);
        th2(i) = thetaf(2);
        th3(i) = thetaf(3);
        th4(i) = thetaf(4);
        th5(i) = thetaf(5);
        th6(i) = [1, t(i)-5*tf, (t(i)-5*tf)^2, (t(i)-5*tf)^3]*tc6; 
    end
        
    xj(i) = (2621*sin(th1(i)))/20000 - (487*cos(th1(i))*cos(th2(i)))/2000 + (921*cos(th5(i))*sin(th1(i)))/10000 + (1707*sin(th2(i) + th3(i) + th4(i))*cos(th1(i)))/20000 - (533*cos(th2(i) + th3(i))*cos(th1(i)))/2500 - (921*cos(th2(i) + th3(i) + th4(i))*cos(th1(i))*sin(th5(i)))/10000;
    yj(i) = (1707*sin(th2(i) + th3(i) + th4(i))*sin(th1(i)))/20000 - (921*cos(th1(i))*cos(th5(i)))/10000 - (487*cos(th2(i))*sin(th1(i)))/2000 - (2621*cos(th1(i)))/20000 - (533*cos(th2(i) + th3(i))*sin(th1(i)))/2500 - (921*cos(th2(i) + th3(i) + th4(i))*sin(th1(i))*sin(th5(i)))/10000;
    zj(i) = 3037/20000 - (533*sin(th2(i) + th3(i)))/2500 - (487*sin(th2(i)))/2000 - (921*sin(th2(i) + th3(i) + th4(i))*sin(th5(i)))/10000 - (1707*cos(th2(i) + th3(i) + th4(i)))/20000;

end

%code for animation
% figure
% view(-38,31)
% for i = 1:length(t)
%     hold on
%     grid on
%     plot3(xj(1:i),yj(1:i),zj(1:i), 'color', 'black', 'LineWidth', 3);
%     axis([-.15,.15,-.15,.15,-.15,.15]);
%     xlabel('xaxis');
%     ylabel('yaxis');
%     zlabel('zaxis');
%     pause(.01)
%     hold off
%     
% end

%plots and figures related to the trajectory planning

figure(1)
grid on
plot3(xj(1:171),yj(1:171),zj(1:171), 'color', 'black', 'LineWidth', 3);
axis([-.15,.15,-.15,.15,-.15,.15]);
xlabel('xaxis');
ylabel('yaxis');
zlabel('zaxis');

figure(3)
grid on
plot(t,xj);
xlabel('t')
ylabel('xaxis')

figure(4)

grid on
plot(t,yj);
xlabel('t')
ylabel('yaxis')

figure(5)
grid on
plot(t,zj);
xlabel('t')
ylabel('zaxis')

figure(2)
subplot(3,2,1)
plot(t,th1)
grid on
xlabel('t')
ylabel('theta_1')
title('theta_1 vs t')

subplot(3,2,2)
plot(t,th2)
grid on
xlabel('t')
ylabel('theta_2')
title('theta_2 vs t')

subplot(3,2,3)
plot(t,th3)
grid on
xlabel('t')
ylabel('theta_3')
title('theta_3 vs t')

subplot(3,2,4)
plot(t,th4)
grid on
xlabel('t')
ylabel('theta_4')
title('theta_4 vs t')

subplot(3,2,5)
plot(t,th5)
grid on
xlabel('t')
ylabel('theta_5')
title('theta_5 vs t')

subplot(3,2,6)
plot(t,th6)
grid on
xlabel('t')
ylabel('theta_6')
title('theta_6 vs t')

%function to calculate jacobian
function J = Jacobian(g,h,i,j,k, ~)
th1 = g;
th2 = h;
th3 = i;
th4 = j;
th5 = k;
J = [(2621*cos(th1))/20000 + (921*cos(th1)*cos(th5))/10000 + (487*cos(th2)*sin(th1))/2000 - (533*sin(th1)*sin(th2)*sin(th3))/2500 + (533*cos(th2)*cos(th3)*sin(th1))/2500 - (1707*cos(th2)*cos(th3)*sin(th1)*sin(th4))/20000 - (1707*cos(th2)*cos(th4)*sin(th1)*sin(th3))/20000 - (1707*cos(th3)*cos(th4)*sin(th1)*sin(th2))/20000 + (1707*sin(th1)*sin(th2)*sin(th3)*sin(th4))/20000 + (921*cos(th2)*cos(th3)*cos(th4)*sin(th1)*sin(th5))/10000 - (921*cos(th2)*sin(th1)*sin(th3)*sin(th4)*sin(th5))/10000 - (921*cos(th3)*sin(th1)*sin(th2)*sin(th4)*sin(th5))/10000 - (921*cos(th4)*sin(th1)*sin(th2)*sin(th3)*sin(th5))/10000,                  (cos(th1)*(1707*cos(th2 + th3 + th4) - 921*cos(th2 + th3 + th4 + th5) + 4264*sin(th2 + th3) + 4870*sin(th2) + 921*cos(th2 + th3 + th4 - th5)))/20000,            (cos(th1)*(1707*cos(th2 + th3 + th4) - 921*cos(th2 + th3 + th4 + th5) + 4264*sin(th2 + th3) + 921*cos(th2 + th3 + th4 - th5)))/20000,     (3*cos(th1)*(569*cos(th2 + th3 + th4) - 307*cos(th2 + th3 + th4 + th5) + 307*cos(th2 + th3 + th4 - th5)))/20000, (921*cos(th1)*cos(th2)*cos(th5)*sin(th3)*sin(th4))/10000 - (921*cos(th1)*cos(th2)*cos(th3)*cos(th4)*cos(th5))/10000 - (921*sin(th1)*sin(th5))/10000 + (921*cos(th1)*cos(th3)*cos(th5)*sin(th2)*sin(th4))/10000 + (921*cos(th1)*cos(th4)*cos(th5)*sin(th2)*sin(th3))/10000,                                                                                                                                                                                                             0;
     (2621*sin(th1))/20000 - (487*cos(th1)*cos(th2))/2000 + (921*cos(th5)*sin(th1))/10000 + (533*cos(th1)*sin(th2)*sin(th3))/2500 - (533*cos(th1)*cos(th2)*cos(th3))/2500 + (1707*cos(th1)*cos(th2)*cos(th3)*sin(th4))/20000 + (1707*cos(th1)*cos(th2)*cos(th4)*sin(th3))/20000 + (1707*cos(th1)*cos(th3)*cos(th4)*sin(th2))/20000 - (1707*cos(th1)*sin(th2)*sin(th3)*sin(th4))/20000 - (921*cos(th1)*cos(th2)*cos(th3)*cos(th4)*sin(th5))/10000 + (921*cos(th1)*cos(th2)*sin(th3)*sin(th4)*sin(th5))/10000 + (921*cos(th1)*cos(th3)*sin(th2)*sin(th4)*sin(th5))/10000 + (921*cos(th1)*cos(th4)*sin(th2)*sin(th3)*sin(th5))/10000,                  (sin(th1)*(1707*cos(th2 + th3 + th4) - 921*cos(th2 + th3 + th4 + th5) + 4264*sin(th2 + th3) + 4870*sin(th2) + 921*cos(th2 + th3 + th4 - th5)))/20000,            (sin(th1)*(1707*cos(th2 + th3 + th4) - 921*cos(th2 + th3 + th4 + th5) + 4264*sin(th2 + th3) + 921*cos(th2 + th3 + th4 - th5)))/20000,     (3*sin(th1)*(569*cos(th2 + th3 + th4) - 307*cos(th2 + th3 + th4 + th5) + 307*cos(th2 + th3 + th4 - th5)))/20000, (921*cos(th1)*sin(th5))/10000 - (921*cos(th2)*cos(th3)*cos(th4)*cos(th5)*sin(th1))/10000 + (921*cos(th2)*cos(th5)*sin(th1)*sin(th3)*sin(th4))/10000 + (921*cos(th3)*cos(th5)*sin(th1)*sin(th2)*sin(th4))/10000 + (921*cos(th4)*cos(th5)*sin(th1)*sin(th2)*sin(th3))/10000,                                                                                                                                                                                                             0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0, (921*sin(th2 + th3 + th4 - th5))/20000 + (1707*sin(th2 + th3 + th4))/20000 - (921*sin(th2 + th3 + th4 + th5))/20000 - (533*cos(th2 + th3))/2500 - (487*cos(th2))/2000, (921*sin(th2 + th3 + th4 - th5))/20000 + (1707*sin(th2 + th3 + th4))/20000 - (921*sin(th2 + th3 + th4 + th5))/20000 - (533*cos(th2 + th3))/2500, (921*sin(th2 + th3 + th4 - th5))/20000 + (1707*sin(th2 + th3 + th4))/20000 - (921*sin(th2 + th3 + th4 + th5))/20000,                                                                                                                                                                                         - (921*sin(th2 + th3 + th4 - th5))/20000 - (921*sin(th2 + th3 + th4 + th5))/20000,                                                                                                                                                                                                             0;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                              sin(th1),                                                                                                                                        sin(th1),                                                                                                            sin(th1),                                                                                                                                                                                                               sin(th2 - th1 + th3 + th4)/2 + sin(th1 + th2 + th3 + th4)/2, cos(th5)*sin(th1) - cos(th1)*cos(th2)*cos(th3)*cos(th4)*sin(th5) + cos(th1)*cos(th2)*sin(th3)*sin(th4)*sin(th5) + cos(th1)*cos(th3)*sin(th2)*sin(th4)*sin(th5) + cos(th1)*cos(th4)*sin(th2)*sin(th3)*sin(th5);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                0,                                                                                                                                                             -cos(th1),                                                                                                                                       -cos(th1),                                                                                                           -cos(th1),                                                                                                                                                                                                               cos(th2 - th1 + th3 + th4)/2 - cos(th1 + th2 + th3 + th4)/2, cos(th2)*sin(th1)*sin(th3)*sin(th4)*sin(th5) - cos(th2)*cos(th3)*cos(th4)*sin(th1)*sin(th5) - cos(th1)*cos(th5) + cos(th3)*sin(th1)*sin(th2)*sin(th4)*sin(th5) + cos(th4)*sin(th1)*sin(th2)*sin(th3)*sin(th5);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                1,                                                                                                                                                                     0,                                                                                                                                               0,                                                                                                                   0,                                                                                                                                                                                                                                                     -cos(th2 + th3 + th4),                                                                                                                                                                                -sin(th2 + th3 + th4)*sin(th5)];


end

%function for inverse kinematics
function thet = ThetaValues(a,b,c,d,e,f,g,h,i,j,k,l)

syms nx ny nz sx sy sz ax ay az Px Pz
sym('Py')
syms th1 th2 th3 th4 th5 th6 real

%DH parameter table
DHTABLE = [0,0,th1,.15185;
           pi/2,0,th2,0;
           0,-.2435,th3,0;
           0,-.2132,th4,.13105;
           (pi/2),0,th5,.08535;
           -(pi/2),0,th6,.0921;];

nx = a;
ny = b;
nz = c;
sx = d;
sy = e;
sz = f;
ax = g;
ay = h;
az = i;
Px = j;
Py = k;
Pz = l;
%transfer matrix    
TDH = [nx, sx, ax, Px;
       ny, sy, ay, Py;
       nz, sz, az, Pz;
       0,  0,  0,  1];


%position of the 5th joint 
Px5 = TDH(1,4) - DHTABLE(6,4)*TDH(1,3);
Py5 = TDH(2,4) - DHTABLE(6,4)*TDH(2,3);

s1 = sqrt(Px5^2 + Py5^2 - DHTABLE(4,4)^2);

%calculation of theta1 value
theta1 = atan2(Py5, Px5) + atan2(DHTABLE(4,4),s1);                     

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
    s2 = sqrt(TDH(1,1)^2+TDH(1,2)-m3^2);
    
    %calculation of theta6 value
    theta6 = atan2(TDH(1,1),TDH(1,2)) + atan2(s2,m3);
    m4 = -(TDH(1,4)*cos(theta1) + TDH(2,4)*sin(theta1) + DHTABLE(6,4)*cos(theta234)*sin(theta5) - DHTABLE(5,4)*sin(theta234));
    m5 = -(TDH(3,4) - DHTABLE(1,4) + DHTABLE(6,4)*sin(theta5) + DHTABLE(5,4)*cos(theta234));
    s3 = (m4^2 + m5^2 - .2132^2 - .2435^2)/(2*.2132*.2435);
    
    %calculation of theta3 value
    theta3 = atan2(sqrt(1-s3^2),s3);
    
    %calculation of theta2 value
    theta2 = atan2(((.2132*cos(theta3)+.2435)*m5 - .2132*sin(theta3)*m4),((.2132*cos(theta3)+.2435)*m4 + .2132*sin(theta3)*m5));
    
    %calculation of theta4 value
    theta4 = theta234 -theta2 - theta3;
    
    thet = [theta1 theta2 theta3 theta4 theta5 theta6];
end
end