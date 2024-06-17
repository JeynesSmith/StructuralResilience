% This script simulates the behaviour for a range of networks based on
% either the governing equations in the original form or in matrix form. In
% this script we perturb the network by adjusting the stimulus, and
% simulating how the network will react. This script is used to generate
% Figures 1,4,6, and Supplementary Figures 2 and 4. 

%% Figure 1
init = 2.*rand(1,4);
S = [1 2 4];
tend = 100;
k = [0.68 0 0.5 0.28 0 0.99 0 0.58 0.69 0 0 0 0 0 0.75];
Plotting_T = 0;
Plotting_U = [init S(1)];

% simulate
for i = 1:3
[t,u] = ode45(@(t,u) odesys_Fig1(t,u,k,S(i)),[(i-1)*tend i*tend],Plotting_U(end,1:4));
Plotting_U = [Plotting_U; u, S(i).*ones(length(t),1)];
Plotting_T = [Plotting_T; t];
end

% plot
figure, hold on
plot(Plotting_T,Plotting_U(:,1),'-','linewidth',2,'Color',[1 0.5 0])
plot(Plotting_T,Plotting_U(:,3),'linewidth',2.5,'Color',[0 0.5 1])
plot(Plotting_T,Plotting_U(:,5),'--','linewidth',2,'Color','k')
legend('Non-Resilient Population (X)','Resilient Population (Z)','Disturbance')
xlabel('Time'),ylabel('Abundance')
set(gca,'FontSize',15)

%% Figure 4a
tend = 100;
S = [1 2 4 8];
load data\Figure4Paras Aa ra da
init = [0.68,2.78,3.89,0];
Plotting_T = 0;
Plotting_U = [init S(1)];

% simulate
for i = 1:4
[t,u] = ode45(@(t,u) odesys_Fig6(t,u,ra,Aa,da,S(i)),[(i-1)*tend i*tend],Plotting_U(end,1:end-1));

Plotting_U = [Plotting_U; u, S(i).*ones(length(t),1)];
Plotting_T = [Plotting_T; t];
end

% plot
figure, hold on
plot(Plotting_T,Plotting_U(:,1),'-','linewidth',2,'Color',[0 0.25 0.25])
plot(Plotting_T,Plotting_U(:,2),'-','linewidth',2,'Color',[0 0.25 0.75])
plot(Plotting_T,Plotting_U(:,3),'linewidth',2.5,'Color',[1 0.5 0])
plot(Plotting_T,Plotting_U(:,5),'--','linewidth',2,'Color','k')
legend('X','Y','Z','Perturbation')
xlabel('Time'),ylabel('Abundance')
set(gca,'FontSize',15)

%% Figure 4b
tend = 100;
S = [1 2 4 8];
load data\Figure4Paras Ab rb db
init = [3.68,0.40,3.40,0];
Plotting_T = 0;
Plotting_U = [init S(1)];

% simulate
for i = 1:4
[t,u] = ode45(@(t,u) odesys_Fig6(t,u,rb,Ab,db,S(i)),[(i-1)*tend i*tend],Plotting_U(end,1:end-1));

Plotting_U = [Plotting_U; u, S(i).*ones(length(t),1)];
Plotting_T = [Plotting_T; t];
end

% plot
figure, hold on
plot(Plotting_T,Plotting_U(:,1),'-','linewidth',2,'Color',[0 0.25 0.25])
plot(Plotting_T,Plotting_U(:,2),'-','linewidth',2,'Color',[0 0.25 0.75])
plot(Plotting_T,Plotting_U(:,3),'linewidth',2.5,'Color',[1 0.5 0])
plot(Plotting_T,Plotting_U(:,5),'--','linewidth',2,'Color','k')
legend('X','Y','Z','Perturbation')
xlabel('Time'),ylabel('Abundance')
set(gca,'FontSize',15)

%% Figure 6
S = [1 2 3 2 1];
init = 2.*rand(1,4);
tend = 100;
figure, hold on
for j = 1:4
load Data\SmallSystem
if j==2
% Network B - need SmallSystem
A = [A zeros(3,1); zeros(1,4)];r=[r;0];d=[d;0];init=[init,init(1)];
A(3,4)=A(3,1);A(3,1)=0; A(4,1)=A(1,2)*1.5; A(4,4)=A(1,1);
elseif j==3
% Network C - need SmallSystem
A(1,3) = 0.2; A(3,3) = -0.9;
elseif j==4
% Network D - need SmallSystem
A(2,1) = 0.2;A(2,2) = -0.5; 
end
Plotting_T = 0;
Plotting_U = [init S(1)];

% simulate
for i = 1:5
[t,u] = ode45(@(t,u) odesys_Fig6(t,u,r,A,d,S(i)),[(i-1)*tend i*tend],Plotting_U(end,1:end-1)); 

Plotting_U = [Plotting_U; u, S(i).*ones(length(t),1)];
Plotting_T = [Plotting_T; t];
end

% plot
plot(Plotting_T,Plotting_U(:,3),'linewidth',2.5)
end
plot(Plotting_T,Plotting_U(:,5),'--','linewidth',2,'Color','k')
legend('Network A','Network B','Network C','Network D','Perturbation')
xlabel('Time'),ylabel('Abundance')
set(gca,'FontSize',15)

%% Supp. Figure 4
init = 2.*rand(1,4);
S = [1 2 3 2 1];
tend = 100;
load Data\ParasMv2_coexist3.mat
Plotting_T = 0;
Plotting_U = [init S(1)];

% simulate
for i = 1:5
[t,u] = ode45(@(t,u) odesys_SFig4(t,u,k,S(i)),[(i-1)*tend i*tend],Plotting_U(end,1:4));

Plotting_U = [Plotting_U; u, S(i).*ones(length(t),1)];
Plotting_T = [Plotting_T; t];
end

% plot
figure, clf, hold on
plot(Plotting_T,Plotting_U(:,4),'-','linewidth',2,'Color',[0 0.25 0.25])
plot(Plotting_T,Plotting_U(:,2),'-','linewidth',2,'Color',[0 0.25 0.75])
plot(Plotting_T,Plotting_U(:,3),'linewidth',2.5,'Color',[1 0.5 0])
plot(Plotting_T,Plotting_U(:,5),'--','linewidth',2,'Color','k')
legend('Y_{2B}','Y_{2A}','Z','Perturbation')
xlabel('Time'),ylabel('Abundance')
set(gca,'FontSize',15)

%% Supp. Figure 2
init = 2.*rand(1,4);
S = [1 2 3 2 1];
tend = 100;
k = [0.5 0.4 0.3 0.6 0.45 0.55 0.35];
Plotting_T = 0;
Plotting_U = [init S(1)];

% simulate
for i = 1:5
[t,u] = ode45(@(t,u) odesys_SFig2(t,u,k,S(i)),[(i-1)*tend i*tend],Plotting_U(end,1:4));

Plotting_U = [Plotting_U; u, S(i).*ones(length(t),1)];
Plotting_T = [Plotting_T; t];
end

% plot
figure, clf, hold on
plot(Plotting_T,Plotting_U(:,1),'-','linewidth',2,'Color',[0 0.25 0.25])
plot(Plotting_T,Plotting_U(:,2),'-','linewidth',2,'Color',[0 0.25 0.75])
plot(Plotting_T,Plotting_U(:,3),'linewidth',2.5,'Color',[1 0.5 0])
plot(Plotting_T,Plotting_U(:,5),'--','linewidth',2,'Color','k')
legend('X','Y','Z','Perturbation')
xlabel('Time'),ylabel('Abundance')
set(gca,'FontSize',15) 



%% functions
function eqns = odesys_Fig1(t,u,k,S)
eqns = zeros(4,1);
% Species I M O N, S is strength of stimulus
% k = [r2 r4 a11 a12 a14 a23 a24 a31 a32 a34 a41 a42 a43 a44 d]
eqns(1) = -k(3)*u(1)^2 + k(4)*u(1)*u(2);
eqns(2) = k(1)*u(2) - k(6)*u(2)*u(3);
eqns(3) = k(9)*u(2)*u(3) - k(8)*u(1)*u(3) - k(15)*S*u(3);
eqns(4) = 0;
end
function eqns = odesys_Fig6(t,u,r,A,d,S)
eqns = zeros(length(r)+1,1);
% Species I M O N, S is strength of stimulus
eqns(1:end-1) = (r + A*u(1:end-1) + d*S).*u(1:end-1);
eqns(end) = 0;
end
function eqns = odesys_SFig4(t,u,k,S)
eqns = zeros(4,1);
% Species I M O N, S is strength of stimulus
% k = [r2 r4 a11 a12 a14 a23 a24 a31 a32 a34 a41 a42 a43 a44 d]
eqns(1) = -k(3)*u(1)^2 + u(1)*(k(4)*u(2) + k(5)*u(4));
eqns(2) = k(1)*u(2) - k(6)*u(2)*u(3) - k(7)*u(2);
eqns(3) = u(3)*(k(9)*u(2)+k(10)*u(4)) - k(8)*u(1)*u(3) - k(15)*S*u(3);
eqns(4) = k(2)*u(4) - k(13)*u(4)*u(3) + k(7)*u(2);
end
function eqns = odesys_SFig2(t,u,k,S)
eqns = zeros(4,1);
% Species I M O, S is strength of stimulus
% k = [a11 a12 d1 r2 a21 a33 a31]
eqns(1) = -k(1)*u(1)^2 + k(2)*u(1)*u(2) - k(3)*S*u(1);
eqns(2) = k(4)*u(2)-k(5)*u(1)*u(2);
eqns(3) = -k(6)*u(3)^2+k(7)*u(1)*u(3);
eqns(4) = 0;
end