% Ce script Matlab décrit une façon d'automatiser la production 
% de résultats, lorsqu'on doit faire une série de simulations 
% en variant un des paramètres d'entrée.

%% Paramètres %%
%%%%%%%%%%%%%%%%

workingfolder = '.\'; % Chemin d'accès au code compilé
binfilename = 'Schroedinger'; % Nom de l'exécutable
ndigit = 8; % Précision utilisée dans les fichiers d'input
nsimul = 1; % Nombre de simulations à faire

Ninters     = repmat(300, 1, nsimul);
tfin        = repmat(1000.0,1,nsimul);
xL          = repmat(-100.0, 1, nsimul);
xR          = repmat(100.0, 1, nsimul);
omega       = repmat(0.02,1,nsimul);
delta       = repmat(0.0,1,nsimul);
%delta       = [0,2,4,6,8,10,12,14,16,18];
x0          = repmat(0.0,1,nsimul);
n           = repmat(16,1,nsimul);
sigma_norm  = repmat(0.05,1,nsimul);
%sigma_norm  = repmat(sqrt(1/(omega(end)))/(xR(end)-xL(end)),1,nsimul);
dt          = repmat(1, 1, nsimul);
%dt          = linspace(0.05,0.5,nsimul);
%dt          = [0.05,0.08,0.10,0.15];

infilename = 'dt'; % Spécifier ici le nom du paramètre à scanner
inputparam = [Ninters;tfin;xL;xR;omega;delta;x0;n;sigma_norm;dt];
inputparam_string = {'Ninters','tfin','xL','xR','omega','delta','x0','n','sigma_norm','dt'};
nparams = length(inputparam_string);

%% Simulations %%
%%%%%%%%%%%%%%%%%

% NB: DO NOT MODIFY THIS (unless you really want to...)

% Loop on the scanned parameter
% For each value, create an input file,
% run the simulation and store in an output file.
% Both input and output file names are stored in
% lists (cell arrays).

fnameinput_list = cell(1, nsimul); % init empty list
fnameoutput_list = cell(1, nsimul); % init empty list

for ii = 1:nsimul
    % loop on scanned parameter
    
    % create the file name
    filename = [infilename, num2str(eval([infilename,'(ii)'])), '.dat'];
    
    % store the input/output file names
    fnameinput_list{ii} = ['inp', filename]; % add the prefix and store
    fnameoutput_list{ii} = ['out', filename]; % add the prefix and store
    
    % create the input data file
    fid = fopen([workingfolder, fnameinput_list{ii} ], 'wt' ); % create or overwrite (empty file, text mode)
    % fill the file
    for jp = 1:nparams
        fprintf(fid, ['%-1s', '=','%.', num2str(ndigit), 'g\n'], inputparam_string{jp}, inputparam(jp,ii));
    end
    fprintf(fid, ['output=./', '%-s\n'], fnameoutput_list{ii});
    fclose(fid);
    display(fnameinput_list{ii})
    display(fnameoutput_list{ii})
    % run the simulation (on Windows OS, do not specify the working folder)
    eval(['!', workingfolder, binfilename, ' ', fnameinput_list{ ii}]);
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des résultats de toutes les simulationsT=[];
if(strcmp(infilename,'Ninters'))
    
elseif(strcmp(infilename,'dt'))
    
end
   
for ii = 1:nsimul
    if(strcmp(infilename,'delta'))
        AnalyseSchroedinger(fnameoutput_list{ii});
    elseif(strcmp(infilename,'dt'))
        [xmoy,pmoy,incert_x,incert_p,E,t,Pgauche,Pdroite,V,x,psi2]=AnalyseSchroedinger(fnameoutput_list{ii},omega);
%         xmoy_end(ii)=xmoy(end-1)+(xmoy(end)-xmoy(end-1))/(t(end)-t(end-1))*(tfin(end)-t(end-1));
%         pmoy_end(ii)=pmoy(end-1)+(pmoy(end)-pmoy(end-1))/(t(end)-t(end-1))*(tfin(end)-t(end-1));
%         incert_x_end(ii)=incert_x(end-1)+(incert_x(end)-incert_x(end-1))/(t(end)-t(end-1))*(tfin(end)-t(end-1));
%         incert_p_end(ii)=incert_p(end-1)+(incert_p(end)-incert_p(end-1))/(t(end)-t(end-1))*(tfin(end)-t(end-1));
%         Pgauche_end(ii)=Pgauche(end-1)+(Pgauche(end)-Pgauche(end-1))/(t(end)-t(end-1))*(tfin(end)-t(end-1));
%         Pdroite_end(ii)=Pdroite(end-1)+(Pdroite(end)-Pdroite(end-1))/(t(end)-t(end-1))*(tfin(end)-t(end-1));
%         x_classique=sqrt(2*E./omega(end)^2).*sin(omega(end)*t);
%         x_classique_end(ii)=x_classique(end-1)+(x_classique(end)-x_classique(end-1))/(t(end)-t(end-1))*(tfin(end)-t(end-1));
%         close all;
    end
end

if(strcmp(infilename,'delta'))
    
elseif(strcmp(infilename,'dt'))
    
end


%% Figure %%

% figure(1)
% plot(dt.^2,xmoy_end,'b.')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<x>(t_{end})')

% figure(2)
% plot(dt.^2,pmoy_end,'b.')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<p>(t_{end})')
% 
% figure(3)
% plot(dt.^2,xmoy_end,'b.')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<\Deltax>(t_{end})')
 
% figure(4)
% plot(dt.^2,incert_p_end,'b.')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<\Deltap>(t_{end})')

figure(5)
plot(t,Pgauche,t,Pdroite,t,Pgauche+Pdroite)
set(gca,'FontSize',16)
grid on
xlabel('t')
ylabel('P')
legend('P_{x<0}','P_{x>0}','P_{tot}','Location','Best')
 
figure(6)
plotyy(t,E,t,Pgauche+Pdroite)
set(gca,'FontSize',16)
grid on
xlabel('t')
ylabel('E')

% figure(7)
% plot(dt.^2,Pgauche_end,'b.')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('p_{x<0}(t_{end})')

% figure(8)
% plot(dt.^2,Pdroite_end,'b.')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('p_{x>0}(t_{end})')

% x_conv=abs(xmoy_end-x_classique_end);

% figure(10)
% plot(dt.^2,x_conv,'b.')
% set(gca,'FontSize',16)
% xlabel('\Deltat^2')
% ylabel({'$|\left<x\right>-x_{Newton}|(t_{end})$'},'Interpreter','latex')

figure (11)
[X,T] = meshgrid(x,t);
pcolor(X,T,psi2)
shading interp
colormap jet
c = colorbar;
xlabel('x [m]')
ylabel('t [s]')
ylabel(c,'|\psi|^2')
