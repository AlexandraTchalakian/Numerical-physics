%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some variables

workingfolder = './'; % Path to the folder that contains the binary 
% of the code and the simulations (must end with /). Here we execute the
% script in the same folder.
binfilename = 'Aiguille'; % Name of the binary executable
ndigit = 8; % Precision used in numerical to string conversion for input data file name
solver='StormerVerlet';

% Define parameters: scanned and constant ones
% NB: scanned parameter must always be the first in these lines.
dt = [0.001,0.003,0.005,0.008,0.1]; %here, scan on dt, other parameters are kept constant.
nsimul = length( dt ); %number of points in the scan.
Ig     = repmat(0.4, 1, nsimul);
mu     = repmat(0.5, 1, nsimul);
B0     = repmat(1.0, 1, nsimul);
B1     = repmat(0.0, 1, nsimul);
nu     = repmat(0.0, 1, nsimul);
nDtParT = repmat(50, 1, nsimul);
tFin   = repmat(40.0,1, nsimul);
omega  = repmat(sqrt(0.5/0.4), 1, nsimul);
theta0 = repmat(0, 1, nsimul);
vtheta0= repmat(1e-6, 1, nsimul);
omega0= sqrt((mu.*B0)./Ig);
Emanal=-mu.*B0.*cos(0.000001)

infilename = { 'dt' }; % Variable name whose value will be included in the simulation file name
inputparam = [Ig; mu; B0; B1; nu; dt; nDtParT; tFin; omega; theta0; vtheta0];
inputparam_string = {'Ig', 'mu', 'B0', 'B1', 'nu', 'dt', 'nDtParT', 'tFin', 'omega', 'theta0', 'vtheta0'};
nparams = length(inputparam_string);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    %Pnc=data(:,5);
    
    soltheta= theta0(ii)*cos(omega0(ii)*t)+vtheta0(ii)*sin(omega0(ii)*t)/omega0(ii);
    solvtheta=-theta0(ii)*omega0(ii)*sin(omega0(ii)*t)+vtheta0(ii)*cos(omega0(ii)*t);
    dtheta(ii)=abs(theta(end)-soltheta(end))/abs(soltheta(end));
    dthetadot(ii)=abs(thetadot(end)-solvtheta(end))/abs(solvtheta(end));
    
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
    
    figure(3)
    %subplot(m,n,2)
    hold on
    plot(t,Emec-Emec(1))
    %set(gca,'yscale','log')
    ylabel('\Delta E_{mec}')
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
    
    dtconv(ii)=t(2)-t(1);
    thetafinalconv(ii)=theta(end);
    dEmec(ii)=abs(Emec(end)-Emanal(end))/abs(Emanal(end));
    
end

save('workspace');

    
figure(1)
hold on
%plot(t,soltheta,'g--')

figure(2)
hold on
plot(t,solvtheta,'g--')

figure(4)
hold on
plot(soltheta,solvtheta,'g--')

ptheta=polyfit(log(dtconv),log(dtheta),1)
vtheta=polyval(ptheta,log(dtconv));
pvtheta=polyfit(log(dtconv),log(dthetadot),1)
vvtheta=polyval(pvtheta,log(dtconv));
pEm=polyfit(log(dtconv),log(dEmec),1)
vEm=polyval(pEm,log(dtconv));

figure
plot(dtconv.^2,thetafinalconv,'x')
xlabel('\Delta t[s]')
ylabel('\theta(t_{fin})')

figure
loglog(dtconv,exp(vtheta),'g-',dtconv,dtheta,'bx')
xlabel('dt [s]')
ylabel('dtheta [rad]')

figure
loglog(dtconv,exp(vvtheta),'g-',dtconv,dthetadot,'bx')
xlabel('dt [s]')
ylabel('d\theta/dt [rad/s]')

figure
plot(dtconv,dEmec,'bx');
xlabel('dt [s]')
ylabel('dEmec [J]')
