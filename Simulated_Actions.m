% This script creates Figure 5, a simulation demonstrating how knowledge of
% perfect resilience can be used to manipulate a network and save a species
% from extinction. In this system, species A alters the network away from a
% perfect resilience network with species O being the target. When the
% network is perturbed, species O will go extinct. By performing an action
% which suppresses species A to extinction, the network match a perfect
% resilience network, and species O is able to recover and persist despite
% furhter perturbations.

% initialise system and parameters
init = 2.*rand(1,4);
load ActionParameters
Plotting_T = [];
Plotting_U = [];
tend = 3;
FinalT = 10*tend;

% simulate initially
S=1; A=0;
[t,u] = ode45(@(t,u) odesys(t,u,k,S,A),[0 tend],init);
Plotting_U = [Plotting_U; u, S.*ones(length(t),1)];
Plotting_T = [Plotting_T; t];

% apply a perturbation
S=3; A=0;
[NoPertT,NoPertU] = ode45(@(t,u) odesys(t,u,k,1,0),[tend FinalT],Plotting_U(end,1:4));
[t,u] = ode45(@(t,u) odesys(t,u,k,S,A),[tend 2*tend],Plotting_U(end,1:4));
Plotting_U = [Plotting_U; u S.*ones(length(t),1)];
Plotting_T = [Plotting_T; t];

% apply action
S=3; A=10;
[NoActT,NoActU] = ode45(@(t,u) odesys(t,u,k,S,0),[2*tend FinalT],Plotting_U(end,1:4));
[t,u] = ode45(@(t,u) odesys(t,u,k,S,A),[2*tend 6*tend],Plotting_U(end,1:4));
Plotting_U = [Plotting_U; u S.*ones(length(t),1)];
Plotting_T = [Plotting_T; t];

% apply perturbation
S=5; A=10;
[t,u] = ode45(@(t,u) odesys(t,u,k,S,A),[6*tend FinalT],Plotting_U(end,1:4));
Plotting_U = [Plotting_U; u S.*ones(length(t),1)];
Plotting_T = [Plotting_T; t];


% plot
Plotting_U=Plotting_U.*40; % scaling
NoPertU=NoPertU.*40; % scaling
NoActU=NoActU.*40; % scaling

figure(2), clf, hold on
plot(Plotting_T,Plotting_U(:,3),'linewidth',2,'Color',[0 0.5 1]) % plot target
plot(Plotting_T,Plotting_U(:,4),'linewidth',2,'Color',[1 0.5 0]) % plot actioned species
plot(Plotting_T,Plotting_U(:,5),'linewidth',2,'Color','k') % plot stimulus
plot(NoActT,NoActU(:,3),'--','LineWidth',2,'Color',[0 0.25 0.5]) % plot no action
plot(Plotting_T,Plotting_U(:,1:2),'.-') % plot other species
legend('Species of concern','Acted Upon Species','Perturbation','No Action','I','M')
xlabel('Time'),ylabel('Abundance')
set(gca,'FontSize',15)


function eqns = odesys(t,u,k,S,A)
eqns = zeros(4,1);
% Species I M O N, S is strength of stimulus, A is action to remove or add
% species
% k = [r2 r4 a11 a12 a14 a23 a24 a31 a32 a34 a41 a42 a43 a44 d]
eqns(1) = -k(3)*u(1)^2 + k(4)*u(1)*u(2);
eqns(2) = k(1)*u(2) - k(6)*u(2)*u(3) - k(7)*u(2)*u(4);
eqns(3) = k(9)*u(2)*u(3) - k(8)*u(1)*u(3) - k(15)*S*u(3) - k(10)*u(3)*u(4);
eqns(4) = -k(13)*u(3)*u(4) - k(14)*u(4)^2 + k(11)*u(1)*u(4) + k(12)*u(2)*u(4) - A*u(4);
end