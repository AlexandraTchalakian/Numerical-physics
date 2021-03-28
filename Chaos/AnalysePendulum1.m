%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
nDtParT = [100,150,200,300,400];
nsimul = length( nDtParT); %number of points in the scan.
dt =repmat(-1,1,nsimul) %here, scan on dt, other parameters are kept constant.
Ig     = repmat(0.4, 1, nsimul);
mu     = repmat(0.5, 1, nsimul);
B0     = repmat(1.0, 1, nsimul);
B1     = repmat(3.0, 1, nsimul);
nu     = repmat(0.08, 1, nsimul);
tFin   = repmat(2000.0,1, nsimul);
omega  = repmat(2.24, 1, nsimul);
%theta0 = repmat(0, 1, nsimul);
%vtheta0= repmat(1e-6, 1, nsimul);
theta0 = repmat(1e-8, 1, nsimul);
vtheta0= repmat(1e-6, 1, nsimul);
omega0= sqrt((mu.*B0)./Ig);
Emanal=-mu.*B0.*cos(0.000001);

infilename = { 'nDtParT' }; % Variable name whose value will be included in the simulation file name
inputparam = [Ig; mu; B0; B1; nu; dt; nDtParT; tFin; omega; theta0; vtheta0];
inputparam_string = {'Ig', 'mu', 'B0', 'B1', 'nu', 'dt', 'nDtParT', 'tFin', 'omega', 'theta0', 'vtheta0'};
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
    %fprintf( fid, [ '%.', num2str( ndigit ), 'g\n' ], inputparam( :, ii ) );
    fprintf(fid,['solver=','%-s\n'], solver);
    for jp = 1 : nparams
        fprintf( fid, [ '%-1s', '=','%.', num2str( ndigit ), 'g\n' ], inputparam_string{jp}, inputparam( jp, ii ) );
    end
    fprintf( fid, ['outputPath=./', '%-s\n'], fnameoutput_list{ ii } );
    fclose( fid );
    fnameinput_list{ ii }
    fnameoutput_list{ ii }
    %run the simulation
    %eval( [ '!cp ', fnameinput_list{ ii}, ' configuration.in' ] );
    %eval( [ '!', workingfolder, binfilename ] );
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
    
    dt(ii)=2*pi/(nDtParT(ii)*omega(ii))
    %dtconv(ii)=t(2)-t(1);
    %thetafinalconv(ii)=theta(end);
    
    
       
    figure(1)
    hold on
    plot(dt(ii),t(end), '+')
    
    
    
    figure(2)
    %m=3;
    %n=1;
    %subplot(m,n,1)
    hold on
    plot(t,theta)
    ylabel('\theta [rad]')
    xlabel('t[s]')
   
    
    figure(3)
    hold on
    plot(t,thetadot)
    ylabel('d\theta/dt [rad/s]')
    xlabel('t[s]')
  
    figure(4)
    %subplot(m,n,3)
    hold on
    plot(mod(theta(1:nDtParT(ii):end),2*pi),thetadot(1:nDtParT(ii):end),'.')
    xlabel('\theta [rad]')
    ylabel('d\theta/dt [rad/s]')
    
    figure(5)
    plot(theta(1:nDtParT(ii):end),thetadot(1:nDtParT(ii):end),'.')
    
        Emecini = [];
    %subplot(m,n,2)
    
    for i=1:length(Emec)
    Emecini(i) = Emec(1);
    end
    
    figure(6)
    hold on
    plot(t,Emec-Emecini.' - Enc)
    ylabel('\Delta E_{mec} - Enc')
    xlabel('t[s]')
    
    E=Emec-Emecini.'-Enc;



thetafinalconv(ii)=abs(theta(end));
Econv(ii)=abs(E(end));

end
figure (7)
plot(dt.^2,Econv,'x')
xlabel('\Delta t^2 [s^2]')
ylabel('E_{final}')

figure (8)
plot(dt.^2,thetafinalconv,'x')
xlabel('(\Delta t^2 [s^2])')
ylabel('\theta_{final}')

    figure(12)
    hold on
    plot(t,abs(fft(theta)))
    plot(t,abs(fft(thetadot)))
    set(gca,'yscale','log')
    ylabel('abs(fft)')
    xlabel('t[s]')

save('workspace');