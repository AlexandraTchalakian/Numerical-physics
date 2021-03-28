% Ce script Matlab décrit une façon d'automatiser la production 
% de résultats, lorsqu'on doit faire une série de simulations 
% en variant un des paramètres d'entrée.

%% Paramètres %%
%%%%%%%%%%%%%%%%

workingfolder = './'; % Chemin d'accès au code compilé
binfilename = 'Schroedinger'; % Nom de l'exécutable
ndigit = 8; % Précision utilisée dans les fichiers d'input
nsimul =1; % Nombre de simulations à faire

% V0 = 2*(0.7854)^2;
     %omega_0 = sqrt(8*20/(200)^2);
% E0 = 0.15226;
n           = repmat(16,1,nsimul)
Ninters     = repmat(300, 1, nsimul);
tfin        = repmat(1000.0,1,nsimul);
xL          = repmat(-100.0, 1, nsimul);
xR          = repmat(100.0, 1, nsimul);
omega       = repmat(0.02,1,nsimul);
delta       = repmat(10,1,nsimul)
x0          = repmat(0,1,nsimul);
sigma_norm  = repmat(sqrt(1/(omega(1))*1/(xR-xL)),1,nsimul);
dt          = repmat(0.1, 1, nsimul);
%dt          = linspace(0.05,1,nsimul);

infilename = 'sigma_norm'; % Spécifier ici le nom du paramètre à scanner
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
   
err_sol=zeros(1,nsimul);

for ii = 1:nsimul
    if(strcmp(infilename,'Ninters'))
        AnalyseSchroedinger(fnameoutput_list{ii});
    elseif(strcmp(infilename,'sigma_norm'))
        [Pdroite,t,V,incert_x]=AnalyseSchroedinger_2(fnameoutput_list{ii});
         %close all;
       K(ii) = max(Pdroite);
       ii
figure (12)
hold on
plot(t,incert_x)
    end
end

if(strcmp(infilename,'Ninters'))
    
elseif(strcmp(infilename,'dt'))
    

end
% for ii=1:length(xmoy)
% xclass(ii) = (sqrt((2*E0))/omega_0)*sin(omega_0*t(ii));
% pclass(ii) = (sqrt((2*E0))/omega_0)*omega_0*cos(omega_0*t(ii));
% end




% figure(1)
% plot(dt.^2,xmoy_end,'bx')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<x>(t_{end})')
% 
% figure(2)
% plot(dt.^2,pmoy_end,'bx')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<p>(t_{end})')
% 
% figure(3)
% plot(dt.^2,xmoy_end,'bx')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<\Deltax>(t_{end})')
% 
% figure(4)
% plot(dt.^2,incert_p_end,'bx')
% set(gca,'FontSize',16)
% xlabel('\Delta t^{2}')
% ylabel('<\Deltap>(t_{end})')
% 
% figure(5)
% plot(t,Pgauche,t,Pdroite,t,Pgauche+Pdroite)
% set(gca,'FontSize',16)
% grid on
% xlabel('t')
% ylabel('P')
% legend('P_{x<0}','P_{x>0}','P_{tot}','Location','Best')
% 
% figure(6)
% plot(t,E)
% set(gca,'FontSize',16)
% grid on
% xlabel('t')
% ylabel('E')
% 
% figure(8)
% plot(t,xmoy)
% set(gca,'FontSize',16)
% grid on
% ylabel('x_{moy}')
% xlabel('t')
% figure(9)
% plot(t,pmoy)
% set(gca,'FontSize',16)
% grid on
% ylabel('p_{moy}')
% xlabel('t')
% 
% figure (10)
% plot(t,xclass,'.')
% hold on
% plot(t,xmoy)
% set(gca,'FontSize',16)
% grid on
% ylabel('x')
% xlabel('t')
% legend ('Newton','Schroedinger')
% 
% 
% figure (11)
% plot(t,pclass)
% hold on;
% plot(t,pmoy)
% set(gca,'FontSize',16)
% grid on
% ylabel('p')
% xlabel('t')
% legend ('Newton','Schroedinger')
% 


%figure (17)
% for ii=1:n
% plot(n,K(ii,1000));
% hold on;
% end
%plot((((2*pi.*n)/(xR(1)-xL(1))).^2)/(V((end+1)*0.5)),K)
% set(gca,'FontSize',16);
% xlabel('E/V0');
% ylabel('Probabilit� passage � droite');
% 
