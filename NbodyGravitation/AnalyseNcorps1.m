format LONGG

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
binfilename = 'Gravitation'; % Name of the binary executable
ndigit = 18; % Precision used in numerical to string conversion for input data file name
% Define parameters: scanned and constant ones
% MODIFY this according to your needs.
% NB: scanned parameter must always be the first in these lines.

 %here, scan on dt, other parameters are kept constant.
 %lamnda = linspace(0.1,3.13,1000)
 %epsilon = logspace(2,3,15);
 %epsilon=linspace(100,7000,15);
 %epsilon=[1000,500,100,10,1,0.1,0.001,0.0001]
 %dt=linspace(50,1000,15);
dt = [900,700,500,200,100,90,60,50,40,30];
nsimul = length(dt); %number of points in the scan.

%dt = repmat(100,1,nsimul)
tfin = repmat(150000,1,nsimul);
G = repmat(6.674e-11,1,nsimul);
Mt = repmat(5.9736e+24,1,nsimul);
Ma = repmat(1e+5,1,nsimul);
Ml = repmat(0,1,nsimul);
Rt = repmat(6378.1e+3,1,nsimul);
Rl = repmat(1738e+3,1,nsimul);
alpha = repmat(1.4*((2*pi)/360),1,nsimul);
epsilon = repmat(100,1,nsimul);


ra0x = repmat(6e+8,1,nsimul);
ra0y = repmat(0,1,nsimul);
rt0x = repmat(0,1,nsimul);
rt0y = repmat(0,1,nsimul);
rl0x = repmat(384748e+3,1,nsimul);
rl0y= repmat(0,1,nsimul);

va0x = repmat(-5e+3*cos(1.4*((2*pi)/360)),1,nsimul);
va0y = repmat(5e+3*sin(1.4*((2*pi)/360)),1,nsimul);

vt0x = repmat(0,1,nsimul);
vt0y = repmat(0,1,nsimul);
vl0x = repmat(0,1,nsimul);
vl0y= repmat(0,1,nsimul);

adapt = repmat(false,1,nsimul);
d = 384748e+3;
v=sqrt(va0x.^2+va0y.^2);
hmint=(-2*G(1)*Mt(1)+sqrt(4*G(1)^2*Mt(1)^2+4*(v(1)^2-2*G(1)*Mt(1)/ra0x(1))*(ra0x(1)*v(1)*sin(alpha(1)))^2))/(2*(v(1)^2-2*G(1)*Mt(1)/ra0x(1)))

infilename = { 'dt' }; % Variable name whose value will be included in the simulation file name
inputparam = [dt; tfin; G; Mt; Ma; Ml; Rt; Rl; alpha;epsilon; ra0x;ra0y;rt0x;rt0y;rl0x;rl0y;va0x;va0y;vt0x;vt0y;vl0x;vl0y;adapt];
inputparam_string = {'dt', 'tfin', 'G', 'Mt', 'Ma', 'Ml', 'Rt', 'Rl', 'alpha','epsilon', 'ra0x','ra0y','rt0x','rt0y','rl0x','rl0y','va0x','va0y','vt0x','vt0y','vl0x','vl0y','adapt'};
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
    %fprintf(fid,['solver=','%-s\n'], solver);
    for jp = 1 : nparams
        fprintf( fid, [ '%-1s', '=','%.', num2str( ndigit ), 'g\n' ], inputparam_string{jp}, inputparam( jp, ii ) );
    end
    fprintf( fid, ['outputPath=./', '%-s\n'], fnameoutput_list{ ii } );
    fclose( fid );
    fnameinput_list{ ii }
    fnameoutput_list{ ii }
    %run the simulation
    %eval( [ '!cp ', fnameinput_list{ ii}, ' configuration.in' ] );%eval -> écrit tout ça dans a commande matlab : !fait écrire dans la commande linux
    %eval( [ '!', workingfolder, binfilename ] );
    % On Linux platforms, uncomment the previous 2 lines and comment the following 2 lines
    % On Windows platforms, comment the previous 2 lines and uncomment the following 2 lines. Execute the present script in the workingfolder directory.
    eval( [ '!copy ', fnameinput_list{ ii}, ' configuration.in' ] );
    eval( [ '!', 'Gravitation' ] );
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MODIFY this according to your needs.

%for each file, load the data and, as an example, plot all the simulations
%on the same plot

n = [];
dEmec =[];
step = [];
dhmin = [];
hmin = [];
for ii = 1 : nsimul
   
    %cs=[];
    filename=[fnameoutput_list{ ii }]
    data = load( filename ); %watch out: {} for lists, [] for vectors
    
  t=data(:,1);
Xt = data(:,2);
Yt = data(:,3);
Xa = data(:,4);
Ya = data(:,5);
Xl = data(:,6);
Yl = data(:,7);
Vxt = data(:,8);
Vyt = data(:,9);
Vxa = data(:,10);
Vya = data(:,11);
Vxl = data(:,12);
Vyl = data(:,13);
dt_2= data(:,14);
h = data(:,15);
Emec=data(:,16);

t2 = linspace(1,10,1001);

n(ii) = length(dt);

hbis=sqrt(Xa.^2+Ya.^2);

[M,index]= min(abs(h))
varier = 3
ind_min = index-varier
ind_max = index+varier
%hmins=[h(index-2),h(index-1),h(index),h(index+1),h(index+2)];
t_precis = linspace(t(ind_min),t(ind_max),100*varier);
h_precis = interp1(t(ind_min:ind_max),h(ind_min:ind_max),t_precis,'spline');
hmin_real(ii) = min(h_precis)
%h_mins=interp1([t(index-2),t(index-1),t(index),t(index+1),t(index+2)],hmins,t_precis,'spline');
% p=polyfit([t(index-1), t(index), t(index+1)],[hmins(1),hmins(2),hmins(3)],2);
% tinter=t(index-1):(t(index+1)-t(index-1))/10000:t(index+1);
% h_min(ii)=min(p(1)*tinter.^2+p(2)*tinter+p(3));
 dhmin(ii)=abs(hmin_real(ii)-hmint);


figure (1)
hold on
p1 = plot(Xa,Ya);
hold on;
p2 = plot(Rt(1)*cos(t2),Rt(1)*sin(t2),'b');
xlabel('X [m]');
ylabel('Y [m]');
legend([p1(1),p2(1)],'Astéroïde','Terre');
axis equal;


%figure (4)
%plot ((1/n).^4,dt);


figure (2)
hold on
plot(t,dt_2,'d');
xlabel('t [s]');
ylabel('\Delta t [s]');

step(ii) = (1/length(t));
dEmec(ii) = abs(Emec(end)-Emec(1));
end

pEmec=polyfit(log(dt),log(dEmec),1)
vEmec=polyval(pEmec,log(dt));
phmin=polyfit(log(dt),log(dhmin),1)
vhmin=polyval(phmin,log(dt));

figure (5)
loglog(dt,exp(vEmec),'g-',dt,dEmec,'bx')
grid on;
xlabel('\Delta t [s]')
ylabel('\Delta E_{mec} [J]')

figure (6)
loglog(dt,exp(vhmin),'g-',dt,dhmin,'bx')
grid on;
xlabel('\Delta t [s]')
ylabel('\Delta h_{min} [m]')

% pEmec=polyfit(log(step),log(dEmec),1)
% vEmec=polyval(pEmec,log(step));
% phmin=polyfit(log(step),log(dhmin),1)
% vhmin=polyval(phmin,log(step));
% 
% figure (7)
% loglog(step,exp(vEmec),'g-',step,dEmec,'bx')
% grid on
% xlabel('\Delta (N_{step}^{-1}) [s]')
% ylabel('\Delta E_{mec} [m]')
% 
% figure (8)
% loglog(step,exp(vhmin),'g-',step,dhmin,'bx')
% grid on
% xlabel('\Delta (N_{step}^{-1}) [s]')
% ylabel('\Delta h_{min} [m]')