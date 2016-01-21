l = 0.5;

file = 'data100lami.csv';
[Ul,   Yl] = reducedata(file,2,6,3,2);

file = 'data100turb0.csv';
[Ut1,  Yt1] = reducedata(file,2,10,2);

file = 'data5600turb0.csv';
[Ut2,  Yt2] = reducedata(file,2,10,2);

%%

figure(1000)
plot([0; Yl ]/l,[0; Ul ],'k-');grid on;hold on;
plot([0; Yt1]/l,[0; Ut1],'k-+');
plot([0; Yt2]/l,[0; Ut2],'k--');

xlabel('1-y/h [-]');
ylabel('u/ub [-]');
legend('laminar','turb re=100','turb re=1000','Location','southeast');

set(gcf,'color','w');