clear all
close all

%Rate of sampling
t = 0:8:1000; %coarse
t = 0:3:1000;  %fine

%Time
T = 0:1:1000;
l = length(t);

%Initial coefficiants. Damping coeff = 0.05
m1 = 10; c1 = 0.77; k1 = 6;

%Choice of decay function
[dm, dk] = decay_fn(t);

%JFF
M = (1+dm).*m1;
%M = ones([1,l]).*m1;
C = c1*ones([1,l]);
%K = (1+dk).*k1;
K = ones([1,l]).*k1;

%Initial Values
w0 = sqrt(K(1)/M(1));
z0 = C(1)/(2*sqrt(K(1)*M(1)));

%Evolution of Values
wn = sqrt(K./M);
zeta=C./(2*sqrt(M.*K));

wd=wn.*sqrt(1-zeta.^2);

wn1 = -1.*wn.*zeta;
wd1 = wd;

rl = (wn1(1) - wn1)./w0;
im = (wd1(1) - wd1)./w0;

%Encoding error:
[rl,im] = error_incorp(rl,im,0,0.025);



%Only Mass Evolutions
delta_m = (im.*(2-im))./((1-im).^2);
[kernel,basis] = optimizer(t,delta_m);
gpMdl1 = fitrgp(t',delta_m,"KernelFunction",kernel,"BasisFunction",basis,OptimizeHyperparameters="auto");
[ypred1,~,yint1] = predict(gpMdl1, T');


%{
%Only Stiffness Evolution
delta_k = -1*(im.*(2-im));
[kernel,basis] = optimizer(t,delta_k);
gpMdl2 = fitrgp(t',delta_k,"KernelFunction",kernel,"BasisFunction",basis,OptimizeHyperparameters="auto");
[ypred2,~,yint2] = predict(gpMdl2, T');
%}

fig = figure(1);
fig.Position(3) = fig.Position(3)*4;
tiledlayout(1,1,'TileSpacing','compact')

nexttile
hold on
plot(T,ypred1,'color','#D3522C', 'LineWidth',1.5);
plot(t,dm,"LineStyle", "--","Color", "Black");
patch([T';flipud(T')],[yint1(:,1);flipud(yint1(:,2))],'k','FaceAlpha',0.05);
hold off
title('\Delta_m')
xlabel('Normalised Time'); ylabel('-\Delta_m (t_s)');
legend({'Prediction','Actual Data','95% Confidence'},'Location','best')

%{
nexttile
hold on
plot(T,ypred2,'color','#D3522C', 'LineWidth',1.5);
plot(t,dk,"LineStyle", "--","Color", "Black");
patch([T';flipud(T')],[yint2(:,1);flipud(yint2(:,2))],'k','FaceAlpha',0.05);
hold off
title('\Delta_k')
xlabel('Normalised Time'); ylabel('-\Delta_k (t_s)');
legend({'Prediction','Actual Data','95% Confidence'},'Location','best')
%}










%Both Mass and Stiffness Evolution
%{
delta_m = -1*(rl./(rl+z0));
[kernel,basis] = optimizer(t,delta_m);
gpMdl1 = fitrgp(t',delta_m,"KernelFunction",kernel,"BasisFunction",basis,OptimizeHyperparameters="auto");
[ypred1,~,yint1] = predict(gpMdl1, T');

delta_k = ((z0*rl.^2) - (1-2*z0^2)*im + im.^2*z0^2)./(z0+rl);
[kernel,basis] = optimizer(t,delta_k);
gpMdl2 = fitrgp(t',delta_k,"KernelFunction",kernel,"BasisFunction",basis,"OptimizeHyperparameters","auto");
[ypred2,~,yint2] = predict(gpMdl2, T');



fig = figure(2);
fig.Position(3) = fig.Position(3)*4;
tiledlayout(2,1,'TileSpacing','compact')

nexttile
hold on
plot(T,ypred1,'color','#D3522C', 'LineWidth',1.5);
plot(t,dm,"LineStyle", "--","Color", "Black");
patch([T';flipud(T')],[yint1(:,1);flipud(yint1(:,2))],'k','FaceAlpha',0.05);
hold off
title('\Delta_m')
xlabel('Normalised Time'); ylabel('-\Delta_m (t_s)');
legend({'Prediction','Actual Data','95% Confidence'},'Location','best')


nexttile
hold on
plot(T,ypred2,'color','#D3522C', 'LineWidth',1.5);
plot(t,dk,"LineStyle", "--","Color", "Black");
patch([T';flipud(T')],[yint2(:,1);flipud(yint2(:,2))],'k','FaceAlpha',0.05);
hold off
title('\Delta_k')
xlabel('Normalised Time'); ylabel('-\Delta_k (t_s)');
legend({'Prediction','Actual Data','95% Confidence'},'Location','best')
%}

