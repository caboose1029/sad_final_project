clear;
clc;
close all

%initial Conditions:
global G phi theta psi Mp R A B C alpha_ beta_ gamma_
%Constants
G = 6.6743e-20;   %km^3kg^-1s^-2
Re = 6378.1; % km

% ind = 1;
% st = 5;
% for i = 0:st:360
%     for j = 0:st:360

%Orientation
phi = deg2rad(10);
theta = deg2rad(0);
psi = deg2rad(0);

%angular velocity
wx = 0;
wy = 10;
wz = 0;

%Masses
Mb = 1000;        %kg
Mp = 5.972e24;    %kg (Earth)

%Orbital Distance
R = [0;0;1000+Re];   %km

%Inertia Tensor
A = 400;
B = 100;
C = 400;


%Rotation matrix from {a} to {c} and {a} to {s}
aDCMc = Rxyz([phi, theta],[2,1]);
aDCMs = Rxyz([phi, theta, psi],[2,1,2]);

%Quaternion from {a} to {c}
aqc = r2q(aDCMc);

%Euler axis and principal angle {a} to {c}
[aUc,PA] = eul(aDCMc);

%Gravity gradient components
alpha_ = A*(1-3*aDCMs(1,3)^2)+B*(1-3*aDCMs(2,3)^2)+C*(1-3*aDCMs(3,3)^2);
beta_ = A*aDCMs(1,1)*aDCMs(1,3)+B*aDCMs(2,1)*aDCMs(2,3)+C*aDCMs(3,1)*aDCMs(3,3);
gamma_ = A*aDCMs(1,2)*aDCMs(1,3)+B*aDCMs(2,2)*aDCMs(2,3)+C*aDCMs(3,2)*aDCMs(3,3);

%Gravity gradient
f_2 = 3/(Mb*norm(R)^2)*[beta_; gamma_; .5*alpha_];

%Force on body due to gravity
F = -G*Mp*Mb/(norm(R)^2)*([0;0;1]+f_2)

%Moment on body due to gravity
M = cross(-R,F)


W = sqrt((G*Mp)/norm(R)^3);


qd=qDot(aqc,[wx,wy,wz],[phi,theta,psi],W)

time = 30;
timesteps = 100;
tspan = linspace(0,time,timesteps);
ic = [wx, wy, wz, aqc'];

tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
[t,x] = ode45(@EOM, tspan, ic,options);

plot(t,x(:,1:3))

%         tz(:,ind) = M;
%         tx(:,ind) = ones(3,1)*i;
%         ty(:,ind) = ones(3,1)*j;
%         tz2((i/st)+1,(j/st)+1) = norm(M);
%         ind = ind+1;

%     end
% end

% surf(reshape(tz(1,:),size(meshgrid(0:st:360))),'FaceAlpha',0.5,'FaceColor','r','EdgeAlpha',0.1)
% hold on;
% surf(reshape(tz(2,:),size(meshgrid(0:st:360))),'FaceAlpha',0.5,'FaceColor','g','EdgeAlpha',0.1)
% surf(reshape(tz(3,:),size(meshgrid(0:st:360))),'FaceAlpha',0.5,'FaceColor','b','EdgeAlpha',0.1)
% legend({'M_x','M_y','M_z'})
% xlabel('phi')
% ylabel('psi')
% zlabel('M')
% 
% figure
% surf(tz2','FaceAlpha',0.5,'EdgeAlpha',0.1)
% xlabel('phi')
% ylabel('psi')
% zlabel('M')

%% Functions

%Equation of motion intitial conditions: ic = [ wx, wy, wz, q1, q2, q3, q4 ]
function dxdt = EOM(t,ic)
global G phi theta psi Mp R A B C alpha_ beta_ gamma_

om_sq = (G*Mp)/norm(R)^3;

phi =   ic(1)*t;
theta = ic(2)*t;
psi =   ic(3)*t;

w1 = ic(1)*cos(psi)+ic(3)*sin(psi);
w2 = ic(2);
w3 = -ic(1)*sin(psi)+ic(3)*cos(psi);

dxdt(1) = (-3*om_sq*gamma_ - (C-B)*ic(2)*ic(3))/A;   % wx dot
dxdt(2) = (3*om_sq*beta_ - (A-C)*ic(3)*ic(1))/B;     % wy dot
dxdt(3) = (-(B-A)*ic(1)*ic(2))/C;                    % wz dot

dxdt(4) =  .5*( .5*(1+cos(phi)+cos(theta)+cos(phi)*cos(theta))^.5 *w1 +( ( sin(phi)*(cos(theta)+1)*(w3+sqrt(om_sq)*sin(theta))-sqrt(om_sq)*cos(theta)*sin(theta)*sin(phi) ) / ( 2*(1+cos(phi)+cos(theta)+cos(phi)*cos(theta))^.5 ) ));                        % q1 dot
dxdt(5) =  .5*( -.5*(1+cos(phi)+cos(theta)+cos(phi)*cos(theta))^.5 *sqrt(om_sq)*cos(theta) +( ( -sin(phi)*sin(theta)*w1-sin(theta)*(cos(phi)+1)*(w3+sqrt(om_sq)*sin(theta)) ) / ( 2*(1+cos(phi)+cos(theta)+cos(phi)*cos(theta))^.5 ) ));                       % q2 dot
dxdt(6) =  .5*( .5*(1+cos(phi)+cos(theta)+cos(phi)*cos(theta))^.5 *(w3+sqrt(om_sq)*sin(theta)) + ( ( -sqrt(om_sq)*cos(theta)*sin(theta)*(cos(phi)+1)-sin(phi)*(cos(theta)+1)*w1 ) / ( 2*(1+cos(phi)+cos(theta)+cos(phi)*cos(theta))^.5 ) ));                        % q3 dot
dxdt(7) = -.5*(  ( w1*sin(theta)*(cos(phi)+1)-sqrt(om_sq)*cos(theta)*sin(phi)*(cos(theta)+1)-(w3+sqrt(om_sq)*sin(theta))*sin(phi)*sin(theta) ) / ( 2*(1+cos(phi)+cos(theta)+cos(phi)*cos(theta))^.5 ) );                        % q4 dot

dxdt=dxdt';
end


%Takes a single or series of rotation angles and axis where 1:x-axis,
%2:y-axis, 3:z-axis. If no axis is specified, defaults to z-axis. If no
%angle is specified, defaults to 0;
function R = Rxyz(ang,ax)

%Ensure vectors are of same length
maxlen = max(length(ang), length(ax));
ang = [ang, zeros(1,maxlen - length(ang))];
ax =  [ax, zeros(1,maxlen - length(ax))];

%Initialize Rotation matrix as Identity matrix
R = eye(3);

%Calculate and apply each rotation
for i = 1:length(ang)
    switch ax(i)
        case 1 %Rotation about x axis
            R_ = [1, 0,           0;
                0, cos(ang(i)), -sin(ang(i));
                0, sin(ang(i)), cos(ang(i))];
        case 2 %Rotation about y axis
            R_ = [cos(ang(i)), 0, sin(ang(i));
                0,            1, 0;
                -sin(ang(i)), 0, cos(ang(i))];
        otherwise %Default rotationa about z axis
            R_ = [cos(ang(i)),-sin(ang(i)),0;
                sin(ang(i)),  cos(ang(i)), 0;
                0,            0,           1];
    end
    
    %Apply rotation to overall rotation matrix
    R = R*R_;
end
end

%Converts rotation matrix to quaternions using provided formulas
function q = r2q(DCM)
q4 = .5*(1+trace(DCM))^(0.5);
q1 = (DCM(3,2)-DCM(2,3))/(4*q4);
q2 = (DCM(1,3)-DCM(3,1))/(4*q4);
q3 = (DCM(2,1)-DCM(1,2))/(4*q4);

q = [q1;q2;q3;q4];
end

%Convers rotation matrix or quaternions to euler axes and principal angle
function [ax,ang] = eul(q)

    %Check if rotation matrix or quaternion
    if size(q) == [3,3]
        q = r2q(q);
    end

    %Calculate principal angle and euler axes
    ang = 2*acos(q(4));
    ax = q(1:3)./sin(ang/2);

end

function q_dot = qDot(q,wb,ang,W)
    w = [wb(1)*cos(ang(3))+wb(3)*sin(ang(3));
        wb(2);
        -wb(1)*sin(ang(3))+wb(3)*cos(ang(3))];


    q_dot4 = (-q(1)*w(1)-q(2)*w(2)-q(3)*(w(3)-wb(3)-W))/2;
    q_dot1 = (q(2)*(w(3)-wb(3)+W)-q(3)*w(2)+q(4)*w(1))/2;
    q_dot2 = (q(3)*w(1)+q(4)*w(2)-q(1)*(w(3)-wb(3)+W))/2;
    q_dot3 = (q(4)*(w(3)-wb(3)-W)+q(1)*w(2)-q(2)*w(1))/2;

    q_dot = [q_dot1;q_dot2;q_dot3;q_dot4];

end

%This is a new comment
