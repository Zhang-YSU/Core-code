%% Dynamic parameter identification
clear;
close all
tic
%% load data
% loadpath = '.\data\data.mat';  %data path
%Torque 1 or current 0
Tau_Cueentflag = 1;
% Symbol analytic output 1 without 0
FLAG_CHART = 0;
%Set the number of friction terms 
nf = 1;  %2 3
if(nf == 0)
    FLAG_FRICTION = 0; % Regardless of friction 0; friction 1
else
     FLAG_FRICTION = 1; 
end
FLAG = 0;%Verify 1 without Verify 0
Selected_Traj = 'you data';
loadpath = ['you data path',Selected_Traj,'name.mat'];
SaveParaPath = ['you data path',Selected_Traj,'_',num2str(nf),'name.mat'];
data = load(loadpath);

q = data.q_f';    
dq = data.dq_f';  
ddq = data.ddq_calculate1_f';  

if Tau_Cueentflag == 0
    efforts = data.current';   
    efforts_f = data.current_f';% efforts_f = data.current_smooth';
else
    efforts = data.force';   
    efforts_f = data.force_f'; % efforts_f = data.force_smooth';
end
%% Calculate the regression matrix corresponding to the kinetic parameters
% M-DH
% alpha=[alpha1;alpha2;alpha3;alpha4;alpha5;alpha6];    
% a=[a1;a2;a3;a4;a5;a6];         
% d=[d1;d2;d3;d4;d5;d6];           
% theta=[theta1;theta2;theta3;theta4;theta5;theta6];
% dh_list = [alpha a d theta];

g = 9.81; 

down_rate = 10;
dt = 0.002*down_rate; 
t = 0:dt:dt*(size(q,2)-1);
n = size(q,2);

% Calculate the linear regression matrix
dof = 6;
if(nf == 0)
    J = zeros(n*dof,6); K = zeros(n*dof,dof*10); Kf=[];
    for i = 1:n
        [J_temp,K_temp,Kf_temp] = Compute_Dynmatrix(q(:,i), dq(:,i), ddq(:,i), g, dh_list, nf);
        for j = 1:dof
            J(n*(j-1)+i,:) = J_temp(j,:); % related to external force
            K(n*(j-1)+i,:) = K_temp(j,:); % related to the inertia parameters of the robot
        end
    end
else
    J = zeros(n*dof,6); K = zeros(n*dof,dof*10); KKf=[];
    for i = 1:n
        [J_temp,K_temp,Kf_temp] = Compute_Dynmatrix(q(:,i), dq(:,i), ddq(:,i), g, dh_list, nf);
        for j = 1:dof
            J(n*(j-1)+i,:) = J_temp(j,:); % related to external force
            K(n*(j-1)+i,:) = K_temp(j,:); % related to the inertia parameters of the robot
            Kf(n*(j-1)+i,:) = Kf_temp(j,:);% related to friction
        end
    end
end
%% Extract the minimum parameter set form
% Eliminate 0 columns ( corresponding to parameters unrelated to the kinetic model )
K_d = K(:,find(sum(abs(K))'>sqrt(eps)));

% QR decomposition is applied to obtain the linear regression matrix K _ min corresponding to the minimum parameter set.
[~,R,P] = qr(K_d);
n_d = rank(K_d);
KK = K_d*P;
K_min = KK(:,1:n_d); 
%% Whether to consider friction
if(FLAG_FRICTION == 1)
    K_min = [K_min,Kf];
end
%% The output kinetic parameter symbolic analytical formula K _ total is the complete parameter set, and K _ min is the minimum parameter set of 36, which is time-consuming.
if(FLAG_CHART == 1)
    K_total = [K,Kf];
    fid=fopen('K_total.txt','w');
    for i = 1:dof_num
        fprintf(fid,['K_total(',num2str(i),',:) = %s;\r'],char(K_total(i,:)));
        fprintf(fid,'####################################\r');
    end
    fclose(fid);
    K_min = [K_min,Kf];
    fid=fopen('K_min.txt','w');
    for i = 1:dof_num
        fprintf(fid,['K_min(',num2str(i),',:) = %s;\r'],char(K_min(i,:)));
        fprintf(fid,'####################################\r');
    end
    fclose(fid);
end
%%
if(FLAG == 1)
    %Verify
    load(VerifyParaPath)
    tau_all = K_min * para_ls;   
    tau_all_m = reshape(tau_all,[n,size(tau_all,1)/n])';
    tau_all_irls = K_min * para_irls;        
    tau_all_irls_m = reshape(tau_all_irls,[n,size(tau_all_irls,1)/n])';
else
    %% LS
    para_ls= (K_min'*K_min)\K_min' * reshape (efforts_f',[],1);
    toc  
    % The accuracy of parameter estimation is verified by backstepping ( that is, the estimated torque and the actual measured torque are verified on the same trajectory ).
    tau_all = K_min * para_ls;   %Using the minimum parameter set obtained by the rigid estimation, the torque is recalculated.
    tau_all_m = reshape(tau_all,[n,size(tau_all,1)/n])';
    %% WLS
    r = var((efforts_f - tau_all_m)');
    sigma_inv = [];
    for i = 1:dof
        sigma_i = ones(n,1)./ r(i);
        sigma_inv = [sigma_inv;sigma_i];
    end
    sigma_inv = diag(sigma_inv);
    para_wls= (K_min'*sigma_inv*K_min)\K_min'*sigma_inv*reshape (efforts_f',[],1); 
    tau_all_wls = K_min * para_wls;        
    tau_all_wls_m = reshape(tau_all_wls,[n,size(tau_all_wls,1)/n])';
    %% IRLS
    para_irls_0=para_ls;
    tau_all_irls = tau_all;
    tau_all_irls_m = tau_all_m;
    while 1
        r = var((efforts_f - tau_all_irls_m)');
        w = [];
        for i = 1:dof
            w_i = ones(n,1)./ r(i);
            w = [w;w_i];
        end
        W = diag(w);  
        para_irls= (K_min'*W*K_min)\K_min'*W*reshape (efforts_f',[],1); 
        tau_all_irls = K_min * para_irls;        
        tau_all_irls_m = reshape(tau_all_irls,[n,size(tau_all_irls,1)/n])';
        error = norm(para_irls-para_irls_0);
        if error<1e-10 
            break;
        end
        para_irls_0=para_irls;
    end
end
%% Error
err_ls = efforts - tau_all_m;  
err_wls = efforts - tau_all_wls_m;
err_irls = efforts - tau_all_irls_m; 
err_f = efforts - efforts_f; %The signal error before and after filtering is assumed to be measurement noise.
%RMSE
RMS_ls = rms(err_ls'); 
RMS_wls = rms(err_wls');
RMS_irls = rms(err_irls'); 
RMS_f = rms(err_f'); 

RMS_compare = [RMS_ls',RMS_wls',RMS_irls',RMS_f'];
%% PLOT
close all
figure
for cnt = 1:6
    subplot(2,3,cnt)
    plot(t,efforts(cnt,:),'g*'),hold on, grid on
    plot(t,efforts_f(cnt,:),'y','linewidth',2.5),hold on, grid on
    plot(t,tau_all_m(cnt,:),'r-','linewidth',1),hold on
    plot(t,tau_all_wls_m(cnt,:),'b-.','linewidth',1),hold on
    plot(t,tau_all_irls_m(cnt,:),'k--','linewidth',1),hold on
    %plot(t,err_ls(cnt,:),'k'),hold on
    if Tau_Cueentflag == 0
        xlabel('t(s)'), ylabel('Current(A)')
        lgd = legend('Current','After filter Current','LS','WLS','IRLS','Location','southeast');
        lgd.Color = 'none';
    else
        xlabel('t(s)'), ylabel('Torque(N.m)')
        lgd = legend('Torque','After filter Torque','LS','WLS','IRLS','Location','southeast');
        lgd.Color = 'none';
    end
    title(['Joint',num2str(cnt)])
end
%err plot
figure
for cnt = 1:6
    subplot(2,3,cnt)
    plot(dq(cnt,:)*180/pi,err_ls(cnt,:),'r*','markersize',1.5),hold on, grid on
    plot(dq(cnt,:)*180/pi,err_wls(cnt,:),'g*','markersize',1.5),hold on, grid on
    plot(dq(cnt,:)*180/pi,err_irls(cnt,:),'b*','markersize',1.5),hold on, grid on
    if Tau_Cueentflag == 0
        xlabel('W(¡ã/s£©'), ylabel('Error(A)')
        lgd = legend('LS','WLS','IRLS','Location','southeast');
    else
        xlabel('W(¡ã/s£©'), ylabel('Error(N.m)')
        lgd = legend('LS','WLS','IRLS','Location','southeast');
    end
    lgd.Color = 'none';
    title(['Before filter Joint',num2str(cnt)])
end
figure
for cnt = 1:6
    subplot(2,3,cnt)
    plot(dq(cnt,:)*180/pi,efforts_f(cnt,:) - tau_all_m(cnt,:),'b*','markersize',1.5),hold on, grid on
    plot(dq(cnt,:)*180/pi,efforts_f(cnt,:) - tau_all_wls_m(cnt,:),'r*','markersize',1.5),hold on, grid on
    plot(dq(cnt,:)*180/pi,efforts_f(cnt,:) - tau_all_irls_m(cnt,:),'g*','markersize',1.5),hold on, grid on
    if Tau_Cueentflag == 0
        xlabel('W(¡ã/s£©'), ylabel('Error(A)')
        lgd = legend('LS','WLS','IRLS','Location','southeast');
    else
        xlabel('W(¡ã/s£©'), ylabel('Error(N.m)')
        lgd = legend('LS','WLS','IRLS','Location','southeast');
    end
    lgd.Color = 'none';
    title(['After filter Joint',num2str(cnt)])
end
%% Save the minimum dynamic parameters and correlation transformation matrix
save(SaveParaPath,'para_irls_0','para_wls','para_ls');