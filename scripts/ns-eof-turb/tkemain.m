[U,   Y] = reducedata(file,colt.BY,colt.ARCL,colt.U1     ,1);
[V,   Y] = reducedata(file,colt.BY,colt.ARCL,colt.U2     ,1);
[NU,  Y] = reducedata(file,colt.BY,colt.ARCL,colt.NUT    ,1);
[EPS, Y] = reducedata(file,colt.BY,colt.ARCL,colt.EPSILON,1);
[TKE, Y] = reducedata(file,colt.BY,colt.ARCL,colt.TKE,1);
[D,   Y] = reducedata(file,colt.BY,colt.ARCL,colt.D,1);

if yplusb
    YY = Y / lplus;
else
    YY = Y;
end

% Dissipation
% plot(YY,-EPS-D,'b');hold on; grid on;

% Production
S12 = ((U(3:end)-U(1:end-2))./(Y(3:end)-Y(1:end-2)));
S12 = S12.*S12.*NU(2:end-1);
plot(YY(2:end-1),S12,sprintf('%s-s',form));hold on; grid on;

% Viscous Diffusion
K12 = (((TKE(3:end)-TKE(1:end-2))./(Y(3:end)-Y(1:end-2)))).*(1/Re+NU(2:end-1));
K12 = (((K12(3:end)-K12(1:end-2))./(Y(4:end-1)-Y(2:end-3))));

plot(YY(3:end-2),K12,sprintf('%s-o',form));

% Turbulent convection
L12 = (TKE(3:end)-TKE(1:end-2))./(Y(3:end)-Y(1:end-2));
L12 = L12.*V(2:end-1);

% plot(YY(2:end-1),-K12,'k');

% Dissipation
plot(YY,-EPS-D,sprintf('%s-*',form));hold on; grid on;
% plot(YY(3:end-2),-K12-S12(2:end-1),sprintf('%s-*',form));

legend('Production','Viscous Diffusion','Dissipation');


set(gcf,'color','w');