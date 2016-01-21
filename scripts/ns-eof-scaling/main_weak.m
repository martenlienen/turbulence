filename = 'channel-weak-intel-128x64-validationevaluation\all.csv';

%%
M = csvread(filename,1,0);

nrs   = M(:,11);

ranks = M(:,2);
total = M(:,12) ./ nrs;
commu = M(:,6)  ./ nrs;
poiss = M(:,9)  ./ nrs;

reduc = total - poiss;

subplot(3,2,1);
% plot(ranks,total(1)*total.^-1);grid on;title('total');
scalingerrorbar(ranks, (total-poiss)./total*100);
title('Share of non PETSc simulation time in total time');
xlabel('Nr. of Processors [-]')
ylabel('Time Proportion [%]');

subplot(3,2,2);
% plot(ranks,commu);grid on;title('Communication');
scalingerrorbar(ranks, commu./(total-poiss)*100);
title('Share of communication in non PETSc simulation time');
xlabel('Nr. of Processors [-]')
ylabel('Time Proportion [%]');

subplot(3,2,3);
% plot(M(:,1),M(:,7));grid on;
scalingerrorbar(ranks, total(1)*total.^-1);
title('Total Speedup');
xlabel('Nr. of Processors [-]')
ylabel('Speedup [-]');

subplot(3,2,5);
% plot(M(:,1),M(:,7));grid on;
scalingerrorbar(ranks, total(1)*total.^-1./ranks*100);
title('Total Efficiency');
xlabel('Nr. of Processors [-]')
ylabel('Efficiency [%]');

subplot(3,2,4);
scalingerrorbar(ranks, reduc(1)*reduc.^-1);
title('Speedup without PETSc');
xlabel('Nr. of Processors [-]')
ylabel('Speedup [-]');

subplot(3,2,6);
scalingerrorbar(ranks, reduc(1)*reduc.^-1./ranks*100);
title('Efficiency without PETSc');
xlabel('Nr. of Processors [-]')
ylabel('Efficiency [%]');

set(gcf,'color','w');

