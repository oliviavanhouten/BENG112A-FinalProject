%% BENG 112A Final Project

clear all;
close all;
clc;

% Compiled by Olivia Van Houten

% Took inspiration from script
% Created by DVJ and Erica Pursell 11/14/2014
% Last modified by DVJ 3/04/2024

%% Loading the data
% Order of data:
% axis_1_load axis_2_load dudx dvdx dudy dvdy

data = load('ProjectBiaxData_v2.txt');      % Choose data set to load

n = length(data(:,1));                      % Length of data

f_1 = data(1:n,1);                          % Load on primary axes
f_2 = data(1:n,2);

dudx = data(1:n,3);                         % Components of displacement gradient
dvdx = data(1:n,4);
dudy = data(1:n,5);
dvdy = data(1:n,6);

freq = 1/200; % 200 Hz is the sampling frequency
tend = (length(data)-1)/200;
t = [0 : freq : tend]; % Setting time vector

%% Deformation gradient F = displacement gradient + I
F = [dudx+ones(n,1) dudy; dvdx dvdy + ones(n,1)]; %% it is a 2 column matrix 
 
b = n/2;
 
F11 = F(1:n,1);     % first upper block 
F12 = F(1:n,2);     % second upper block
F21 = F(n+1:end,1); % first bottom block
F22 = F(n+1:end,2); % second bottom block

% Deformation gradient plots

figure(1);
plot(t,F11,'.-r',t,F22,'.-b',t,F12,'.-c',t,F21,'.-m','MarkerSize',12);
set(gca,'Fontsize',16);
xlabel('time[sec]');
ylabel('Deformation gradient');
title('Deformation Gradient vs Time')
grid on;
legend('F(11)', 'F(22)','F(12)', 'F(21)');

figure(2);
plot(t,F11,'.-r',t,F22,'.-b','MarkerSize',12);
set(gca,'Fontsize',16);
xlabel('time[sec]');
ylabel('Deformation gradient');
title('Normal Strains');
grid on;
legend('F(11)', 'F(22)');

figure(3)
plot(t,F12,'.-r',t,F21,'.-b','MarkerSize',12);
set(gca,'Fontsize',16);
xlabel('time [sec]');
ylabel('Deformation gradient');
title('Shear Strains');
grid on;
legend('F(12)', 'F(21)');

%% Finding S - Nominal - 1st Piola Kirchoff Stress

for i = 1:n
    F_mat{i} = [F11(i) F12(i); F21(i) F22(i)]; % converting F into an array
end

for i = 1:n
    lam{i} = eig(F_mat{i}); % eigenvalues of F
    J{i} = lam{i}(1)*lam{i}(2); % J = det(F)
end

a1 = 2.2301 * 9.89;
a2 = 2.2301 * 9.6;

for i = 1:n             % assume no shear and constant area
    T{i}(1,1) = f_1(i)/a1;
    T{i}(1,2) = 0;
    T{i}(2,1) = 0;
    T{i}(2,2) = f_2(i)/a2;
end

for i = 1:n
    F_inv{i} = pinv(F_mat{i});
end

for i = 1:n
    S{i} = J{i}.*F_inv{i}.*T{i};
end

for i = 1:n
    ST{i} = S{i}';
    ST11(i) = ST{i}(1,1);
    ST12(i) = ST{i}(1,2);
    ST21(i) = ST{i}(2,1);
    ST22(i) = ST{i}(2,2);
end

figure(4)
plot(t,ST11,'.-r',t,ST22,'-b','MarkerSize',12);
set(gca,'Fontsize',20);
xlabel('time [sec]');
ylabel('1st Piola Kirchhoff Stress');
title('Normal Stresses (1PK)');
grid on;
legend('S^T(11)', 'S^T(22)');

%% Finding P - Lagrangian 2nd Piola Kirchoff Stress

for i = 1:n
    P{i} = S{i}.*F_inv{i};
    P11(i) = P{i}(1,1);
    P12(i) = P{i}(1,2);
    P21(i) = P{i}(2,1);
    P22(i) = P{i}(2,2);
end

figure(5)
plot(t,P11,'.-r',t,P22,'-b','MarkerSize',12);
set(gca,'Fontsize',20);
xlabel('time [sec]');
ylabel('2nd Piola Kirchhoff Stress');
title('Normal Stresses (2PK)');
grid on;
legend('P(11)', 'P(22)');

%% Finding E - Green Almonsi Strain 

for i = 1:n
    FT{i} = F_mat{i}';
    C{i} = FT{i} * F_mat{i};
    E{i} = (1/2) * (C{i}-eye(2));
    E11(i) = E{i}(1,1);
    E12(i) = E{i}(1,2);
    E21(i) = E{i}(2,1);
    E22(i) = E{i}(2,2);
end

figure(6)
plot(t,E11,'.-r',t,E22,'-b',t,E12,'-c',t,E21,'-m','MarkerSize',12);
set(gca,'Fontsize',20);
xlabel('time [sec]');
ylabel('Green-Almonsi Strain');
title('Normal and Shear Strains');
grid on;
legend('E(11)', 'E(22)','E(12)','E(21)');

figure(7)
plot(t,E11,'.-r',t,E22,'-b','MarkerSize',12);
set(gca,'Fontsize',20);
xlabel('time [sec]');
ylabel('Green-Almonsi Strain');
title('Only Normal Strains');
grid on;
legend('E(11)','E(22)');

%% Stress vs Strain Graphs

% figure(8)
% plot(E11,P11,'.-r',E22,P22,'-b','MarkerSize',12);
% set(gca,'Fontsize',20);
% xlabel('Strain');
% ylabel('Stress');
% title('Stress vs Strain');
% grid on;
% legend('x-direction stress & strain', 'y-direction stress & strain');

figure(8)
plot(E11,ST11,'.-r',E22,ST22,'-b','MarkerSize',12);
set(gca,'Fontsize',20);
xlabel('Strain');
ylabel('Stress');
title('Stress vs Strain');
grid on;
legend('x-direction stress & strain', 'y-direction stress & strain');

%% Decomposition of F = R U

for i = 1:n
    tf{i} = isdiag(C{i});
    if tf{i} == 1
        U{i} = sqrt(C{i});
    else
        [N{i},C_tilde{i}] = eig(C{i});
        N_T{i} = N{i}';
        U{i} = N{i} * sqrt(C_tilde{i}) * N_T{i};
    end
end

for i = 1:n
    U11(i) = U{i}(1,1);
    U12(i) = U{i}(1,2);
    U21(i) = U{i}(2,1);
    U22(i) = U{i}(2,2);
end

for i = 1:n
    U_inv{i} = pinv(U{i});
    R{i} = F_mat{i} * U_inv{i};
    R11(i) = R{i}(1,1);
    R12(i) = R{i}(1,2);
    R21(i) = R{i}(2,1);
    R22(i) = R{i}(2,2);
end

figure(9)
plot(t,R11,'.-r',t,R22,'-b','MarkerSize',12);
set(gca,'Fontsize',20);
xlabel('time [sec]');
ylabel('Rotation Components');
title('Rotation Tensor');
grid on;
legend('R(11)', 'R(22)');

figure(10)
plot(t,U11,'.-r',t,U22,'-b','MarkerSize',12);
set(gca,'Fontsize',20);
xlabel('time [sec]');
ylabel('Stretch Components');
title('Stretch Tensor');
grid on;
legend('U(11)', 'U(22)');



