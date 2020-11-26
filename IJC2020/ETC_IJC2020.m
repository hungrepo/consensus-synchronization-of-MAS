%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This code is used to generate data in the paper entitled:

%----"Consensus/synchronization of networked nonlinear multiple agent systems with event-triggered communications mechanism 
%----for multi agent nonlinear system", accepted for publication in International Journal of Control, 2020 
%
%----Author: Nguyen Tuan Hung, Antonio Pascoal, Institute System and Robotic, IST, Lisbon
%----Contact: nguyen.hung@tecnico.ulisboa.pt 

%% Important NOTE: Run CVX_SETUP first.
%%%%%%
function ETC_TAC
clear all;
% close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                          Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation time and sampling time
      T=10; Ts=1e-3; 
      number_int=T/Ts;
%% Agents' dynamics
      A = [ 0      1      0      0;
           -2      -1   2   0;
           0       0      0      2;
           1.95    0     -1.95   0];
      B = [0; 1;  0;  0];
      % f = [0; 0   ;  0; -0.333sin(x3)] => Lipchizt constant gamma=0.333 %
      gamma=0.333;
      n=length(A);                  % size of the state
      m=length(B(2,:));             % size of the input
      rng('default');
      x1=-5*rand(4,1);  
      x2= 5*rand(4,1);
      x3= 5*rand(4,1);
      x4= 5*rand(4,1);
      x5=5*rand(4,1);
      x6=5*rand(4,1);
      x7=5*rand(4,1);
      x8=-2*rand(4,1);
      %  x2=x1+0.1*rand(4,1); x3=x1; x4=x1; x5=x1; x6=x1; x7=x1; x8=x1;
      
%% Initialize network
% Adjancy matrix
      Adja = [0  0  0  0  0  1 ;    
              1  0  0  0  0  1 ;   
              0  1  0  0  1  0 ;
              0  0  1  0  0  0 ;
              0  0  1  1  0  0 ;
              0  1  0  0  1  0 ]; 
      N=length(Adja(:,1)); 
      one=ones(N,1);                                            % define vector one
      D=diag(Adja*one);                                         % Degree matrix    
      L=D-Adja;                                                 % Laplacian matrix
%       [U,L_eig] = eig(L);                                            
%       tem = sort(L_eig*one);
      L_eig=eig((L+L')/2);
      tem=sort(L_eig);
      lamda2 = tem(2);
      lamdaN = tem(N);
      one=ones(N,1);
      I_n=eye(n);
      I_N=eye(N);
      W=kron(eye(N)-one*one'/N, I_n);        % projection matrix
      R=eye(N);
%% solve LMI to find P and tau with CVX       
% cvx_begin sdp
% variable P(n,n) symmetric;
% P >= .1*eye(4);
% variable tau; 
% tau > 0;
% variable mu;
% mu  > .1;
% % tau < 5;
% [A*P + P*A'-2*tau*(B*B') +  mu*gamma*I_n   P;...   
%  P                                   -mu*(1/gamma)*I_n] < -10*eye(8);
% cvx_end

% If you dont have cvx use mu tau and P as follows
mu = 95.4558;
tau = 1.1257e+03;
P = [   109.9174 -118.0489   75.1311  -64.9442;
       -118.0489  514.5727   40.7710   37.4653;
        75.1311   40.7710  113.9943  -37.1640;
       -64.9442   37.4653  -37.1640   66.8825] ;
P_inv=inv(P);
   
%% Compute feedback control gains        
        K=-B'*P_inv;
        c=tau/lamda2; 

%% Event trigger threashod function
%  h_i=c0+c1*exp(-c2*t(end));               % exponential form with
       c0=1e-3; c1=5; c2=1;                                                    
%% Compute H_1,H_2,F_1,F_2
RP=kron(R,P_inv);
RP_eig=eig(RP);
tem2 = sort(RP_eig);
H1=A*P + P*A' -2*tau*B*(B)'+mu*gamma*I_n+gamma*P^2/mu;
% H1=A*P + P*A' -2*tau*B*B'+ gamma^2*I_n+ P^2;
H2=-kron(I_N,P)*kron(R,H1)*kron(I_N,P);
H2_eig=eig(H2);
tem3=sort(H2_eig);
F1=RP*W*kron(Adja,-c*B*K);
F2=RP*W*kron(L,-c*B*K);
num1=2*norm(F1,2)*tem2(N*n)*sqrt(N)*c0;
num2=2*norm(F2,2)*tem2(N*n)*sqrt(N)*c0;
dem=tem3(1)*tem2(1);
r1=num1/dem;
r2=num2/dem;
%% Compute lower bound for inter event time tau_1
Term1=B*B'*P_inv;
G=kron(c*L,Term1);
T=kron(Adja,Term1);

cu=c1+c0;
xi=W*[x1;x2;x3;x4;x5;x6];
V0=xi'*RP*xi;
beta_bound=V0/tem2(1);

num3=2*norm(F1,2)*tem2(N*n)*sqrt(N)*cu;
gamma_bound=num3/dem;

u_bar=norm(G)*(beta_bound+gamma_bound)+norm(T)*sqrt(N)*cu;

tau1=log(1+c0*(norm(A)+gamma)/(norm(B)*u_bar))/(norm(A)+gamma);
%% Store data during simulation 
x=cell(N,1);                                                                % For state of N agent
x{1,1}=x1;  x{2,1}=x2;  x{3,1}=x3;  x{4,1}=x4;  
x{5,1}=x5;  x{6,1}=x6; x{7,1}=x7; x{8,1}=x8; 
u=cell(N,1);                                                                % For control input of N agent
Com_Sig=cell(N,1);                                                          % For communication signals of N agent;
e=cell(N,1);                                                                % For estimation errors of N agents;
h=cell(N,1);                                                                % For threashod functions of N agents;
delta=cell(N,1);                                                            % For triggering function of N agents; 
Com_mode= 2;                                                                % 1: continous, 2: ETC with time, 3: ETC with state 
t=[];
V=[];
xi_norm=[];
x_hat=cell(N,1);
% u_hat=cell(N,1);                                                            % Not neccessary
%     for i=1:N 
%         x_hat{i,1}=zeros(n,1);
%         u_hat{i,1}=zeros(m,1);
%     end
x_hat=x;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%                                                      Start  simulation 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
for step=1:number_int
 t(end+1)=(step-1)*Ts;                                                      % Compute currrent time
 x1=x{1,1}(:,end);  x2=x{2,1}(:,end);   x3=x{3,1}(:,end);   
 x4=x{4,1}(:,end);  x5=x{5,1}(:,end);  x6=x{6,1}(:,end);   
 xi=W*[x1;x2;x3;x4;x5;x6];
 V(:,end+1)=xi'*(kron(R,P_inv))*xi;                                         % Compute Lyapunov function
 xi_norm(:,end+1)=sqrt(xi'*xi);
 
%% Check mode of communication mode and update control law
    if Com_mode==1                                                          % continous communication 
       for i=1:N
           tem=0;
           for j=1:N
               tem=tem+c*K*Adja(i,j)*(x{i,1}(:,end)-x{j,1}(:,end));
           end
           u{i,1}(:,end+1)=tem;
       end
    end
%% Time dependent trigger
  if Com_mode==2     
  % Update control input for each agent
        for i=1:N
           tem=0;
           for j=1:N
               % if use control law (4)
               tem=tem+c*K*Adja(i,j)*(x{i,1}(:,end)-x_hat{j,1}(:,end));     % For control law (4) in the paper
               % if use control law (42) 
%              tem=tem+c*K*Adja(i,j)*(x_hat{i,1}(:,end)-x_hat{j,1}(:,end));
           end
           u{i,1}(:,end+1)=tem;
        end
   % Compute estimation error and triggering functions
        for i=1:N
            e{i,1}(:,end+1)=x{i,1}(:,end)-x_hat{i,1}(:,end);                % Compute e_i
            h{i,1}(:,end+1)=c0+c1*exp(-c2*t(end));                          % Compute h_i
            delta{i,1}(:,end+1)=norm(e{i,1}(:,end))-h{i,1}(:,end);          % Compute delta_i
%            delta{i,1}(:,end+1)=norm(c*K*e{i,1}(:,end))-h{i,1}(:,end); 
        end
  % Broadcast signal
        for i=1:N
            if  delta{i,1}(:,end)>= 0
                x_hat{i,1}(:,end)=x{i,1}(:,end);                            
                u_hat{i,1}(:,end+1)=u{i,1}(:,end);
                Com_Sig{i,1}(:,end+1)=1;                                    % Means i broadcasts signal
            else
                Com_Sig{i,1}(:,end+1)=0;                                    % no communication
%                 u_hat{i,1}(:,end+1)=u_hat{i,1}(:,end);
            end
        end     
  end
%% Update agents' state
for k=1:N
%       if (t(end)>100)&&(k==1)
%                x{k,1}(:,end+1)=update_x(x{k,1}(:,end),u{k,1}(:,end)+10*exp(-(t(end)-10)),Ts,t(end));
%       else
              x{k,1}(:,end+1)=update_x(x{k,1}(:,end),u{k,1}(:,end),Ts,t(end));
%      end
%% Update estimated (predicted) agents' state
    x_hat{k,1}(:,end+1)=update_xhat(x_hat{k,1}(:,end),0*u_hat{k,1}(:,end),Ts,t(end));
end
%% 
% eta(:,end+1)=W*x(:,end);
% V(:,end+1)=0.5*eta(:,end)'*eta(:,end);
end
    save_to_base(1);
end
function x_next=update_x(x_current,input,Ts,t)
%     input_x.nSteps = 4;
%     input_x.Ts=Ts;
%     input_x.u=input;
%     input_x.t=t;
%     input_x.x = x_current;
%     output_x= RK4_integrator(@xdot, input_x);
%     x_next = output_x.value;
     [time y]=ode45(@(t,y) xdot(t,y,input), [0, Ts], x_current);
     x_next=y(end,:)';
end
function x_next=update_xhat(x_current,input,Ts,t)
%     input_xhat.nSteps = 4;
%     input_xhat.Ts=Ts;
%     input_xhat.u=input;
%     input_xhat.t=t;
%     input_xhat.x = x_current;
%     output_xhat= RK4_integrator(@xhatdot, input_xhat);
%     x_next = output_xhat.value;
    [time y]=ode45(@(t,y) xhatdot(t,y,input), [0, Ts], x_current);
     x_next=y(end,:)';
end
function dx = xdot(t,x,u)
% model of each agent   
n=length(x); 
      A = [ 0      1      0      0;
           -2      -1     2      0;
           0       0      0      2;
           1.95    0     -1.95   0];
      B = [0; 1;  0;  0];
      f = [0; 0   ;  0; -0.333*sin(x(3))];
      dx=A*x+f+B*u;
end
function dx = xhatdot(t,x,u)
% Model to predict neighbor agents 
    n=length(x); 
      A = [ 0      1      0      0;
           -2      -1   2   0;
           0       0      0      2;
           1.95    0     -1.95   0];
      B = [0; 1;  0;  0];
      f = [0; 0   ;  0; -0.333*sin(x(3))];
      dx=A*x+f+B*u;
end
