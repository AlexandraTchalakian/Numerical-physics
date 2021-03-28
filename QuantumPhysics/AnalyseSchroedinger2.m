function [Pdroite,t,V,incert_x]=AnalyseSchroedinger(fichier)

%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fichier = 'output';
data = load([fichier,'_obs.dat']);
t = data(:,1);
Pgauche = data(:,2);
Pdroite = data(:,3);
E = data(:,4);
xmoy = data(:,5);
x2moy = data(:,6);
pmoy = data(:,7);
p2moy = data(:,8);
incert_x=data(:,9);
incert_p=data(:,10);
data = load([fichier,'_pot.dat']);
x = data(:,1);
V = data(:,2);
psi2 = load([fichier,'_psi2.dat']);

xmoy_end=xmoy(end);
pmoy_end=pmoy(end);
incert_x_end=incert_x(end);
incert_p_end=incert_p(end);
%% Figures %%
%%%%%%%%%%%%%
figure('Name',['Analyse de ' fichier])
subplot(2,2,1)
plot(x,V), hold on
plot(x([1,end]),E(1)*ones(1,2),'--')
grid on
xlabel('x')
ylabel('V')
legend('V(x)','E','Location','Best')

subplot(2,2,2)
[X,T] = meshgrid(x,t);
pcolor(X,T,psi2)
shading interp
colormap jet
c = colorbar;
xlabel('x [m]')
ylabel('t [s]')
ylabel(c,'|\psi|^2')

subplot(2,2,3)
plot(t,Pgauche,t,Pdroite,t,Pgauche+Pdroite)
grid on
xlabel('t')
ylabel('P')
legend('P_{x<0}','P_{x>0}','P_{tot}','Location','Best')

subplot(2,2,4)
plot(t,E)
grid on
xlabel('t')
ylabel('E')

% figure 
% [X,T] = meshgrid(x,t);
% pcolor(X,T,psi2)
% shading interp
% colormap jet
% c = colorbar;
% xlabel('x [m]')
% ylabel('t [s]')
% ylabel(c,'|\psi|^2')
% 
% figure 
% plot(t,Pgauche,t,Pdroite,t,Pgauche+Pdroite)
% grid on
% xlabel('t')
% ylabel('P')
% legend('P_{x<0}','P_{x>0}','P_{tot}','Location','Best')
% 
% figure
% plot(x,V), hold on
% plot(x([1,end]),E(1)*ones(1,2),'--')
% grid on
% xlabel('x')
% ylabel('V')
% legend('V(x)','E','Location','Best')

save('data')