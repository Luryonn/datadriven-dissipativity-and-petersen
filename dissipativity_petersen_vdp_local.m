% % % Data-driven Control - Dissipativity-based with Petersen's Lemma
% % % Author: João Gabriel Napoleão Silva - NEACON - UFC
% % % System: Van der Pol oscillator

clc; clear; close all

%----------------------------- DATA COLLECT ------------------------------%
global delta
delta = 1e-2;   % Disturbance amplitude

x0 = [-0.2; 0.1]; % initial condition
opt = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'NormControl','off');
[t, xsol] = ode45(@vander, [0 3], x0, opt); % simulation time = [0,1]
X00 = xsol';
[X11, ~] = gradient(xsol', 1);

U00 = [];

for i = 1:length(t)
    U00(:,i) = sin(2*pi*t(i)); % control signal
end

Ts = 30; % number of samples

X0 = X00(:,1:Ts);    % selecting first T samples of the experiment
X1 = X11(:,1:Ts);
U0 = U00(:,1:Ts);
t = t(1:Ts);

% Data-driven Dissipativity Polynomial Control
pvar x1 x2
vars = [x1;x2];
x = [x1;x2];

% e = 1;
% f = [x2;
%     -x1 + e*(1-x1^2)*x2];
% 
% g = [0 1];
% 
h = monomials(x,1:2);

m = 1; % u columns
n = 2; % x columns
p = size(h,1); % h length

DELTA = sqrt(Ts*delta)*eye(n); %Disturbance upper bound

% Initialization of regressors
Z = [monomials(x,1:2)];
W = monomials(x,0);

N = length(Z);
M = length(W);

% plot of generated data
figure
subplot(2,1,1)
plot(t,U0,'LineWidth',1.5)
title('Input')
ylabel('u')
grid on

subplot(2,1,2)
plot(t,X0(1,:),'--',t,X0(2,:),'-.','LineWidth',1.5)
legend('x1','x2','Location','best')
title('Open-loop response')
ylabel('x')
xlabel('t')
grid on

% construction of the data matrices Ad, Bd, Cd
Z0 = [X0(1,:); X0(2,:); X0(1,:).^2; X0(1,:).*X0(2,:); X0(2,:).^2];
W0 = U0;

Ad = [Z0;W0]*[Z0;W0]';
Bd = -[Z0;W0]*X1';
Cd = X1*X1'-DELTA*DELTA';

% construction of the data matrices zeta and Q
zeta = -Ad\Bd;
Qd = Bd'*(Ad\Bd)-Cd;

zeta_N = zeta(1:N,:);
zeta_M = zeta((N+1):(N+M),:);
Ad_i = Ad^(-1/2);
Ad_i_N = Ad_i(:,1:N);
Ad_i_M = Ad_i(:,(N+1):(N+M));

% Inputs (design parameters)
beta_V = 1e-6;
beta_T = 1e-6;
n_V = 1;
n_T = 2;
epsi = 1e-8;
alpha = 1e-6*(x1^2+x2^2)*eye(n+m+1);

% Maximum iteration number
k = 1;
kmax = 50;

tic
% STEP 1: Determine V0, Tx, L0, Q0, S0, R0, lambda0
prog = sosprogram(vars);

[prog, Q] = sospolymatrixvar(prog,monomials(x,0),[p,p],'symmetric');
[prog, S] = sospolymatrixvar(prog,monomials(x,0),[p,m]);
[prog, R] = sospolymatrixvar(prog,monomials(x,0),[m,m],'symmetric');
[prog, V] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
[prog, T] = sospolyvar(prog,[monomials(x,2:8)],'wscoeff');
[prog, lambda] = sospolyvar(prog,[monomials(x,0:4)],'wscoeff');
[prog, L] = sospolyvar(prog,[monomials(x,0:4)],'wscoeff');

% SOS constraints
prog = sosineq(prog,V-beta_V*(x1^2+x2^2)^n_V); % V radially unbounded constraint
prog = sosineq(prog,T-beta_T*(x1^2+x2^2)^n_T); % Tx radially unbounded constraint

gradV = [diff(V,x1); diff(V,x2)]; % gradient of Lyapunov function

Sigma1 = [gradV'*zeta_N'*Z+T-h'*Q'*h (1/2)*gradV'*zeta_M'*W-h'*S
        W'*zeta_M*(1/2)*gradV-S'*h -R];
Sigma2 = [(1/2)*Z'*Ad_i_N'; (1/2)*W'*Ad_i_M']*[(1/2)*Ad_i_N*Z (1/2)*Ad_i_M*W];
B_aux = [gradV'*Qd^(1/2);
       zeros(m,n)];

stability = [Sigma1+lambda*Sigma2 B_aux          % data-driven dissipativity-based condition for local stability
               B_aux' -lambda*eye(n)]+alpha*(L-1);
prog = sosineq(prog,-stability);

prog = sosineq(prog,R);     % dissipativity's R constraint

prog = sosineq(prog,lambda-epsi); % petersen's lemma's lambda constraint

% STEP 1 solution
sol = sossolve(prog);

Q0 = double(sosgetsol(sol,Q));
R0 = double(sosgetsol(sol,R));
S0 = double(sosgetsol(sol,S));
L = sosgetsol(sol, L);

T = sosgetsol(sol,T);
delta = S0*(R0\(S0'))-Q0;
min(eig(delta))

if min(eig(delta)) >= 0     % Stop criteria: Delta_c test
    K = -R0\(S0')           % Control gains
    V = sosgetsol(sol,V)    % Lyapunov function 
else
    % STEP 2: Iterative method for V, Tx, alpha, Q, S, R, lambda
    while k <= kmax
        prog = sosprogram(vars);

        [prog, Q] = sospolymatrixvar(prog,monomials(x,0),[p,p],'symmetric');
        [prog, S] = sospolymatrixvar(prog,monomials(x,0),[p,m]);
        [prog, R] = sospolymatrixvar(prog,monomials(x,0),[m,m],'symmetric');
        [prog, V] = sospolyvar(prog,[monomials(x,2:4)],'wscoeff');
        [prog, T] = sospolyvar(prog,[monomials(x,2:8)],'wscoeff');
        [prog, lambda] = sospolyvar(prog,[monomials(x,0:4)],'wscoeff');
        [prog, alpha] = sospolymatrixvar(prog,[monomials(x,0:4)],[n+m+1,n+m+1]);
        
        % SOS constraints
        prog = sosineq(prog,V-beta_V*(x1^2+x2^2)^n_V); % V radially unbounded constraint
        prog = sosineq(prog,T-beta_T*(x1^2+x2^2)^n_T); % Tx radially unbounded costraint

        gradV = [diff(V,x1); diff(V,x2)]; % gradient of Lyapunov function
        
        Sigma1 = [gradV'*zeta_N'*Z+T-h'*Q'*h -h'*S+(1/2)*gradV'*zeta_M'*W
                -S'*h+W'*zeta_M*(1/2)*gradV -R];
        Sigma2 = [(1/2)*Z'*Ad_i_N'; (1/2)*W'*Ad_i_M']*[(1/2)*Ad_i_N*Z (1/2)*Ad_i_M*W];
        B_aux = [gradV'*Qd^(1/2)
               zeros(m,n)];

        stability = [Sigma1+lambda*Sigma2 B_aux          % data-driven dissipativity-based condition for local stability
                        B_aux' -lambda*eye(n)]+alpha*(1-L);
        prog = sosineq(prog,-stability);

        prog = sosineq(prog,R);             % dissipativity's R constraint
        
        prog = sosineq(prog,R0-R);          % increasing delta constraints
        delta_increasing = S*(R0\(S0')) + (R0\S0')'*S' - 2*S0*(R0\(S0')) + Q0 - Q; % - 0.01*eye(p);
        prog = sosineq(prog,delta_increasing);

        prog = sosineq(prog,lambda-epsi); % petersen's lemma's lambda constraint
              
        % STEP 2 solution
        sol = sossolve(prog);

        Q = double(sosgetsol(sol,Q))
        R = double(sosgetsol(sol,R))
        S = double(sosgetsol(sol,S))
        
        T = sosgetsol(sol,T);
        delta = S*(R\(S'))-Q;
        min(eig(delta))
        if min(eig(delta))>= 0 | k == kmax     % Stop criteria: Delta_c test
            K = -R\(S')                        % Control gains
            V = sosgetsol(sol,V)               % Lyapunov function 
            break
        end
        k=k+1
        Q0 = Q; R0 = R; S0 = S;
    end
end
toc

% Phase diagram - Open loop
syms x1 x2
a=5;
[x1,x2] = meshgrid(-a:0.1:a,-a:0.1:a);
e=1;
dx1 = x2;
dx2 = -x1 + e.*(1-x1.^2).*x2;

figure
OL = streamslice(x1,x2,dx1,dx2,1);
set(OL,'LineWidth',1.5);
title('Phase diagram - open-loop')
xlabel('x1'), ylabel('x2')
grid on

% Phase diagram - Closed loop
syms x1 x2
a=6;
[x1,x2] = meshgrid(-a:0.1:a,-a:0.1:a);

Kx = K(1).*x1 + K(2).*x2 + ...
     K(3).*x1.^2 + K(4).*x1.*x2 + K(5).*x2.^2 ;
dx1 = x2;
dx2 = -x1 + e.*(1-x1.^2).*x2 + Kx;

figure
CL = streamslice(x1,x2,dx1,dx2,1);
set(CL,'LineWidth',1.5);
title('Phase diagram - close-loop')
xlabel('x1'), ylabel('x2')
hold on
plot(0,0,'o','LineWidth',3)
grid on

%------------------------------ EXPERIMENT -------------------------------%
function dx = vander(t,x)
    e=1;
    Delta = 1e-3;
    u = sin(2*pi*t);
    dx = [x(2) + sqrt(Delta)*cos(2*pi*0.4*t);
         -x(1) + e*(1-x(1)^2)*x(2) + u + sqrt(Delta)*sin(2*pi*0.4*t)];
end