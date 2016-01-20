clear all
global file Re l form colt

%% Half height of the channel
l = 0.5;

%% Settings
colt = reducedatacols('keturb-de'); colt.BY = colt.H;
file = 'data_2_13750_fine.csv';
Re   = 13750;       tauwb = true; div = 1;
form = 'b'; 
tauwmain();

%% Modify stress plot
figure(2000)
set(gcf,'color','w');

subplot(1,3,1);ylim([0,1]);plot([0 1], [1 0],'k--');
subplot(1,3,2);ylim([0,1]);plot([0 1], [1 0],'k--');

%% Modify w+-y+-plot
figure(3000)

y1 = logspace(0,1.1,10);
semilogx(y1,y1,'k--');

y2 = logspace(0.9,3,10);
semilogx(y2,log(y2)/0.41+5.2,'k--');
set(gcf,'color','w');
xlabel('y+ [-]');ylabel('w+ [-]');

legend(...
    'measurement Re=1000',...
    'measurement Re=5600',...
    'measurement Re=13750',...
    'w+=y+',...
    'w+=log(y+)/0.41+5.2','Location','southeast');