

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scan on parameters from Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ce script Matlab decrit une facon d'automatiser la production 
% de resultats, lorsqu'on doit faire une serie de simulations 
% en variant un des parametres d'entree.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some variables

% MODIFY this according to your file system setup
workingfolder = './'; % Path to the folder that contains the binary 
% of the code and the simulations (must end with /). Here we execute the
% script in the same folder.
binfilename = 'Aiguille'; % Name of the binary executable
ndigit = 8; % Precision used in numerical to string conversion for input data file name
solver='StormerVerlet';

% Define parameters: scanned and constant ones
% MODIFY this according to your needs.
% NB: scanned parameter must always be the first in these lines.

 %here, scan on dt, other parameters are kept constant.
 %0.006,0.01875,0.0375,0.075,0.01875,0.0375,0.075,
 %dt = [0.01,0.03,0.06,0.08,0.1]
 mu = 0.5;
 B0 = 1.0;
 Ig = 0.4;
 C= sqrt((mu*B0)/Ig)
%theta0 =  linspace(-pi+0.1,pi-0.1,100);
%theta0 = [pi/8,pi/4,pi/2,pi/2+pi/10,pi/2+pi/5,pi/2+pi/3]
omega = linspace(C-0.05,C+0.05,20);
nsimul = length(omega); %number of points in the scan.
dt     = repmat(0.01, 1, nsimul);
Ig     = repmat(0.4, 1, nsimul);
mu     = repmat(0.5, 1, nsimul);
B0     = repmat(1.0, 1, nsimul);
B1     = repmat(3.0, 1, nsimul);
nu     = repmat(0.08, 1, nsimul);
nDtParT = repmat(50, 1, nsimul);
tFin   = repmat(600,1, nsimul);
%tFin   = repmat(40.0,1, nsimul);
theta0 = repmat(0,1,nsimul);
vtheta0= repmat(1e-6, 1, nsimul);
omega0= sqrt((mu.*B0)./Ig);
%omega = repmat(0, 1, nsimul);
Emanal=-mu.*B0.*cos(0.000001);

infilename = { 'omega' }; % Variable name whose value will be included in the simulation file name
inputparam = [Ig; mu; B0; B1; nu; dt; nDtParT; tFin; omega;theta0; vtheta0];
inputparam_string = {'Ig', 'mu', 'B0', 'B1', 'nu', 'dt', 'nDtParT', 'tFin', 'omega','theta0', 'vtheta0'};
nparams = length(inputparam_string);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NB: DO NOT MODIFY THIS (unless you really want to...)

%Loop on the scanned parameter
%For each value, create an input file,
%run the simulation and store in an output file.
%Both input and output file names are stored in
%lists (cell arrays).

fnameinput_list = cell( 1, nsimul ); %init empty list
fnameoutput_list = cell( 1, nsimul ); %init empty list

for ii = 1 : nsimul
    
    %loop on scanned parameter
    
    %create the file name
    filename = ''; %init
    for jj = 1 : length( infilename )
        %loop on infilename. Build the string with _NameValue
        filename = [ filename, '_', infilename{ jj }, num2str( eval( [ infilename{ jj }, '( ii )' ] ) ) ]; 
    end
    filename = [ filename, '.dat' ]; %add suffix
    
    %store the input/output file names
    fnameinput_list{ ii } = [ 'inp', filename ]; %add the prefix and store
    fnameoutput_list{ ii }  = [ 'out', filename ]; %add the prefix and store
    
    %create the input data file
    fid = fopen( [ workingfolder, fnameinput_list{ ii } ], 'wt' ); %create or overwrite (empty file, text mode)
    %fill the file
    fprintf( fid, [ '%.', num2str( ndigit ), 'g\n' ], inputparam( :, ii ) );
    fprintf(fid,['solver=','%-s\n'], solver);
    for jp = 1 : nparams
        fprintf( fid, [ '%-1s', '=','%.', num2str( ndigit ), 'g\n' ], inputparam_string{jp}, inputparam( jp, ii ) );
    end
    fprintf( fid, ['outputPath=./', '%-s\n'], fnameoutput_list{ ii } );
    fclose( fid );
    fnameinput_list{ ii }
    fnameoutput_list{ ii }
    %run the simulation
    eval( [ '!cp ', fnameinput_list{ ii}, ' configuration.in' ] );
    eval( [ '!', workingfolder, binfilename ] );
    % On Linux platforms, uncomment the previous 2 lines and comment the following 2 lines
    % On Windows platforms, comment the previous 2 lines and uncomment the following 2 lines. Execute the present script in the workingfolder directory.
    eval( [ '!copy ', fnameinput_list{ ii}, ' configuration.in' ] );
    eval( [ '!', binfilename ] );
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MODIFY this according to your needs.

%for each file, load the data and, as an example, plot all the simulations
%on the same plot
dtnum=[];
Tnum=[];
for ii = 1 : nsimul
    dtnum = [];
    %cs=[];
    filename=[fnameoutput_list{ ii }]
    data = load( filename ); %watch out: {} for lists, [] for vectors
    
    t=data(:,1);
    theta=data(:,2);
    thetadot=data(:,3);
    Emec=data(:,4);
    Pnc=data(:,5);
    Enc = data(:,7);
    %DWnc perdu pour chaque dt : pour obtnir l'intégralité des Wnc perdu il faut faire :
    %sommedt Wnc(dt)
    Enctot = [];
    Enctot(1) = Enc(1);
    
    %for i=2:length(Wnc)
        %Wnctot(i) = Wnctot(i-1)+Wnc(i);
    %end
    
    %boucle permettant de stocker dans un vecteur les valeurs de points
    %séparés d'un periode
    
    
    
    for j=2: length(t)
        if((theta(j-1)<0) )
            if((theta(j)>=0))
            dtnum(end+1) = j;
            end
        end
    end
      
       Tnumsum = dtnum(end)-dtnum(1);
    
   
        P = (Tnumsum)/(length(dtnum)-1) %période moyenne calculée pour une simulation
    
       Tnum(ii) = P;     
         

   
    figure(1)
    hold on
    plot(t,mod(theta2*,pi))
    ylabel('\theta [rad]')
    xlabel('t[s]')
   
    
    figure(1)
    hold on
    plot(t,theta)
    ylabel('\theta [rad]')
    xlabel('t[s]')
    
    figure(2)
    hold on
    plot(t,thetadot)
    ylabel('d\theta/dt [rad/s]')
    xlabel('t[s]')
    
    Emecini = [];
    %subplot(m,n,2)
    
    for i=1:length(Emec)
    Emecini(i) = Emec(1);
    end
    
    figure(3)
    hold on
    plot(t,Emec-Emecini.' - Enc)
    ylabel('\Delta E_{mec} - Enc')
    xlabel('t[s]')
    
    figure(4)
    %subplot(m,n,3)
    hold on
    plot(theta,thetadot)
    xlabel('\theta [rad]')
    ylabel('d\theta/dt [rad/s]')
    
    
    figure(5)
    hold on
    plot(dt(ii),t(end), '+')
    
    figure (9)
    hold on
    xlabel('\Deltat [s]')
    ylabel('Enc [J]')
    plot(t,Enc.')
    
    figure (10)
    hold on
    xlabel ('\Deltat [s]')
    ylabel('Emec')
    plot(t,Emec)
     
   
    dtconv(ii)=t(2)-t(1);
    thetafinalconv(ii)=abs(theta(end));
    thetavfinalconv(ii)=abs(thetadot(end));
    Encfinalconv(ii) = abs(Enc(end));
    Difffinalconv(ii) = abs(Emec(end) - Emecini(end) - Enc(end));
    dEmec(ii)=abs(Emec(end)-Emanal(ii))/abs(Emanal(ii));
    
    
    
end

%transformation de la période en radian pour tous les éléments

%for i=1 : length(Tnum)
    %Tnum(i) = (2*pi)/(Tnum(i));
%end
save('workspace');


figure (6)
plot(dtconv.^2,Wncfinalconv,'x')
xlabel('\Delta t[s]');
ylabel('Wnc [J]');

figure (7)
plot(dtconv.^2,Difffinalconv,'x')
xlabel('\Delta t[s]');
ylabel('\DeltaE_{mec} - Wnc [J]');


figure (8) 
plot(theta0,Tnum,'bx');
xlabel('theta0 [rad]')
ylabel('T [ms]')


