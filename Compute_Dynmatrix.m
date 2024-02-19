function [J,K,Kf] = Compute_Dynmatrix(thetalist, dthetalist, ddthetalist, g, dh_list, nf)
% The dynamic model is linearized and the linear regression matrix is calculated.
% Input : thetalist: nx1,joint angle; dthetalist: nx1,joint angular velocity; dthetalist: nx1,angular acceleration of joints; 
%         g: 1x1,acceleration of gravity; dh_list: nx4, M-DH; nf: 1x1, Number of friction terms
% output: J: nx6, Linear regression matrix related to external force
%         K: nx(10*n), Linear regression matrix related to robot body parameters
%         Kf: nx(nf*n), Linear regression matrix related to friction.
%         nf=1: M1;  nf=2: M2; nf=3: M3; 
%%
n = size(dh_list,1);

alpha = dh_list(:,1);
a = dh_list(:,2);
d = dh_list(:,3);
theta = dh_list(:,4);

Z=[0;0;1];
%Establishment of transformation matrix

theta = theta + thetalist;
T=zeros(4,4,n);R=zeros(3,3,n);P=zeros(3,n);
for i=1:n
    T(:,:,i)=[cos(theta(i))                   -sin(theta(i))                0                  a(i)
              sin(theta(i))*cos(alpha(i))  cos(theta(i))*cos(alpha(i))  -sin(alpha(i))   -d(i)*sin(alpha(i))
              sin(theta(i))*sin(alpha(i))  cos(theta(i))*sin(alpha(i))  cos(alpha(i))    d(i)*cos(alpha(i))
              0                              0                              0                  1];
    R(:,:,i)=T(1:3,1:3,i);
    P(:,i)=T(1:3,4,i);
end

%Forward recursion of kinematics
w0 = zeros(3,1); dw0 = zeros(3,1);
dv0 = [0;0;g];
w = zeros(3,n); dw = zeros(3,n);
dv = zeros(3,n);

%i = 0
w(:,1) = R(:,:,1)' * w0 + dthetalist(1) * Z;
dw(:,1) = R(:,:,1)' * dw0 + cross(R(:,:,1)' * w0, dthetalist(1) * Z) + ddthetalist(1) * Z;
dv(:,1) = R(:,:,1)' * (cross(dw0,P(:,1)) + cross(w0,cross(w0, P(:,1))) + dv0);
for i = 1:n-1
   w(:,i+1) = R(:,:,i+1)' * w(:,i) + dthetalist(i+1) * Z ;
   dw(:,i+1) = R(:,:,i+1)' * dw(:,i) + cross(R(:,:,i+1)' * w(:,i), dthetalist(i+1) * Z)+ ddthetalist(i+1) * Z;
   dv(:,i+1) = R(:,:,i+1)' * (cross(dw(:,i), P(:,i+1)) + cross(w(:,i), cross(w(:,i), P(:,i+1))) + dv(:,i));
end

A = Compute_Amatrix(R,P);
B = Compute_Bmatrix(w,dw,dv);

% J
V = []; 
for i = 1:n
    V_temp = eye(6,6);
    for j = i:n
    V_temp = V_temp * A(:,:,j);
    end
    V = [V;V_temp];
end

J = [];
for i = 1:n
    J = [J;V(6*i,:)];
end

% K
U = [];
for i =1:n-1
    RR = eye(6,6);
    D = [];
    for j = i:n-1
        RR = RR * A(:,:,j);
        D = [D,RR * B(:,:,j+1)];
    end
    U_temp = [zeros(6,(i-1)*10),B(:,:,i),D];
    U = [U;U_temp];
end
U = [U;zeros(6,(n-1)*10),B(:,:,n)]; 

K = [];
for i =1:n
    K = [K;U(6*i,:)];
end
%%
Kf = zeros(n,nf*n);
epsilon = 0.001;
for i = 1:n
    if(nf == 1) %M1
        Kf(i,:) = [zeros(1,nf*(i-1)) [tanh(dthetalist(i)/epsilon) dthetalist(i) ddthetalist(i)] zeros(1,nf*(n-i))];
    elseif(nf == 2) %M2
        Kf(i,:) = [zeros(1,nf*(i-1)) [tanh(dthetalist(i)/epsilon)*(tanh(dthetalist(i)/epsilon)+1)/2  dthetalist(i)*(tanh(dthetalist(i)/epsilon)+1)/2 tanh(dthetalist(i)/epsilon)*(1-tanh(dthetalist(i)/epsilon))/2  dthetalist(i)*(1-tanh(dthetalist(i)/epsilon))/2 ] zeros(1,nf*(n-i))];
    elseif(nf == 3) %M3
        Kf(i,:) = [zeros(1,nf*(i-1)) [tanh(dthetalist(i)/epsilon)*(tanh(dthetalist(i)/epsilon)+1)/2  dthetalist(i)*(tanh(dthetalist(i)/epsilon)+1)/2 tanh(dthetalist(i)/epsilon)*(1-tanh(dthetalist(i)/epsilon))/2  dthetalist(i)*(1-tanh(dthetalist(i)/epsilon))/2 ddthetalist(i)] zeros(1,nf*(n-i))];
    end
end
