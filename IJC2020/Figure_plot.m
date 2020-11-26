close all;
fontsize=12;
COL1= [0 0 0]/255;             % black
COL2= [0 0 255]/255;           % Blue 
COL3= [255 0 0]/255;         % Red
COL4= [0,128,128]/255;           % red 
COL5= [255 0 255]/255;         
COL6= [255 0 20]/255;         
Col=[COL1;COL2;COL3;COL4;COL5;COL6];
f1=figure(1);
f2=figure(2);
f3=figure(3);
f4=figure(4);
n=length(t);
for i=1:N
set(0, 'currentfigure', f1);  %# for figuresp
hold on; grid on;         limit=[-0 10 -inf inf]; axis(limit);
ylabel(' $x_{i,1}$','Interpreter','latex','FontSize', fontsize);
xlabel('t[s]','Interpreter','latex','FontSize', fontsize);
plot(t,x{i,1}(1,1:n),'LineWidth',1);

set(0, 'currentfigure', f2);  %# for figuresp
hold on;grid on;         limit=[-0 10 -inf inf]; axis(limit);
ylabel('$x_{i,2}$','Interpreter','latex','FontSize', fontsize);
xlabel('t[s]','Interpreter','latex','FontSize', fontsize);
plot(t,x{i,1}(2,1:n),'LineWidth',1);

set(0, 'currentfigure', f3);  %# for figuresp
hold on;grid on;         limit=[-0 10 -inf inf]; axis(limit);
ylabel(' $x_{i,3}$','Interpreter','latex','FontSize', fontsize);
xlabel('t[s]','Interpreter','latex','FontSize', fontsize);
plot(t,x{i,1}(3,1:n),'LineWidth',1);

set(0, 'currentfigure', f4);  %# for figuresp
hold on;grid on;         limit=[-0 10 -inf inf]; axis(limit);
ylabel(' $x_{i,4}$','Interpreter','latex','FontSize', fontsize);
xlabel('t[s]','Interpreter','latex','FontSize', fontsize);
plot(t,x{i,1}(4,1:n),'LineWidth',1);
end
set(f1,'position',[0 0 560 200]);
set(f2,'position',[0 0 560 200]);
set(f3,'position',[0 0 560 200]);
set(f4,'position',[0 0 560 200]);


% tmin=t(1); tmax=t(end)+1;
% limit=[tmin tmax -inf inf];
% axis(limit);
% set(gca,'YTick',[0 0.4 0.7 1]);
% set(gca,'Yticklabel',[]); 
 
% legend('Agent 1','Agent 2','Agent 3', 'Agent 4', 'Agent 5');
% Plot communication 
  fig5=figure(5);
  for i=1:1:N
  set(fig5,'position',[0 0 560 220]);
  idx=find(Com_Sig{i,1}(1,1:n)==1);
  plot(t(idx),i*Com_Sig{i,1}(1,idx),'.','MarkerSize',8,'Color',Col(i,:)); hold on;
  end
%   idx=find(Com(2,:)==1);
%   plot(t(idx),1*Com(2,idx),'+','Color','k'); hold on;
%   idx=find(Com(3,:)==1);
%   plot(t(idx),0.8*Com(3,idx),'+','Color','g'); hold on;
%   idx=find(Com(4,:)==1);
%   plot(t(idx),0.6*Com(4,idx),'+','Color','m'); hold on;
%   idx=find(Com(5,:)==1);
%   plot(t(idx),0.4*Com(5,idx),'+','Color','c'); hold on;
%   idx=find(Com(6,:)==1);
%   plot(t(idx),0.2*Com(6,idx),'+','Color','r'); hold on;

   
%   tmin=t(1); tmax=t(end)+1;
   limit=[-0 10 -0 7];
   axis(limit);
%   set(gca,'YTick',[0 0.4 0.7 1]);
  set(gca,'Yticklabel',[]); 
%   grid on
 % title('Broadcast time instants of the agents','Interpreter','latex');
  xlabel('t[s]','Interpreter','latex','FontSize', fontsize);
  ylabel('$t_{i,k}$','Interpreter','latex','FontSize', fontsize);

%  legend('Agent 1','Agent 2','Agent 3', 'Agent 4', 'Agent 5');
 
%% Plot estimation error 
  fig6=figure(6);
  set(fig6,'position',[0 0 560 600]);
  for i=1:N
  subplot(N,1,N-i+1);
  plot(t,h{i,1}(1,1:n),'LineWidth',1); hold on;
  plot(t,h{i,1}(1,1:n)+delta{i,1}(1,1:n),'LineWidth',1,'Color',Col(i,:)); 
  limit=[-0 10 -0 5];    axis(limit);
  lg=legend('$h_1(t)$','$||{\bf e}_2(t)||$' );
  set(lg,'Interpreter','latex');

  end
%   xlabel('t[s]','Interpreter','latex');

%% Compute minimum inter-event time
%   n=10000;
%   MIE=zeros(N,1);
%   for i=1:1:N
%   idx=find(Com_Sig{i,1}(1,1:n)==1);
%   MIE(i)=idx(2)-idx(1);
%   for j=2:length(idx)-2
%       new_interval=idx(j+2)-idx(j+1);
%       if  new_interval< MIE(i)
%       MIE(i)=new_interval;
%       end
%   end
%   end
%   MIE=MIE*Ts;

%% Plot the asymtotic bound of synchronizatoin error 
% fig6=figure(6);
% set(fig6,'position',[0 0 560 250]);
% plot(t,xi_norm,'LineWidth',2);
% hold on;
% plot(t,r1*ones(1,length(t)));
% plot(t,r2*ones(1,length(t)));

