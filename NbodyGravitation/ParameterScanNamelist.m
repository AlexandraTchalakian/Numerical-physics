%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some variables
workingfolder = './'; % Path to the folder that contains the binary 
% of the code and the simulations (must end with /). Here we execute the
% script in the same folder.
binfilename = 'Gravitation'; % Name of the binary executable
ndigit = 18; % Precision used in numerical to string conversion for input data file name
% Define parameters: scanned and constant ones.
% NB: scanned parameter must always be the first in these lines.

 %here, scan on dt, other parameters are kept constant.
 %lamnda = linspace(0.1,3.13,1000)
  %epsilon = linspace(10,500,30);
  %dt=linspace(100,1000,5)
dt = [100];
nsimul = length(dt); %number of points in the scan.

%dt = repmat(1,1,nsimul)
tfin = repmat(300000,1,nsimul);
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

%va0y = repmat(0,1,nsimul);
vt0x = repmat(0,1,nsimul);
vt0y = repmat(0,1,nsimul);
vl0x = repmat(0,1,nsimul);
vl0y = repmat(0,1,nsimul);

adapt = repmat(true,1,nsimul);
d = 384748e+3;


%Exercice 3
rt0x = repmat((-Ml(1)*d)/(Ml(1)+Mt(1)),1,nsimul);
rt0y = repmat(0,1,nsimul);
rl0x = repmat((Mt(1)*d)/(Ml(1)+Mt(1)),1,nsimul); 
rl0y = repmat(0,1,nsimul);

%vt0x = repmat(-Ml(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*sin(1.4*((2*pi)/360)),1,nsimul);
%vt0y = repmat(Ml(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*cos(1.4*((2*pi)/360)),1,nsimul);
%vl0x = repmat(Mt(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*sin(1.4*((2*pi)/360)),1,nsimul);
%vl0y = repmat(-Mt(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*cos(1.4*((2*pi)/360)),1,nsimul);

%vt0x = repmat ();

infilename = { 'dt' }; % Variable name whose value will be included in the simulation file name
inputparam = [dt; tfin; G; Mt; Ma; Ml; Rt; Rl; alpha;epsilon; ra0x;ra0y;rt0x;rt0y;rl0x;rl0y;va0x;va0y;vt0x;vt0y;vl0x;vl0y;adapt];
inputparam_string = {'dt', 'tfin', 'G', 'Mt', 'Ma', 'Ml', 'Rt', 'Rl', 'alpha','epsilon', 'ra0x','ra0y','rt0x','rt0y','rl0x','rl0y','va0x','va0y','vt0x','vt0y','vl0x','vl0y','adapt'};
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


%for each file, load the data and, as an example, plot all the simulations
%on the same plot

n = []

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
dt = data(:,14);
hmin = data(:,15);
%Emec=data(:,15);


t2 = linspace(1,10,1000);
n(ii) = length(dt);

figure (1)
p1 = plot(Xa,Ya,'r');
%p2 = plot(Xt,Yt,'b');
hold on;
%p1 = plot(Xl,Yl,'r');
%hold on;
p2 = plot(Rt(1)*cos(t2),Rt(1)*sin(t2),'b');
axis equal;
xlabel('X [m]');
ylabel('Y [m]');
legend([p1(1),p2(1)],'asteroine','Terre');


figure (2)
plot (Xl,Yl,'b');

figure (3)
plot (Vxt,Vyt,'b');
hold on;
plot (Vxl,Vyl,'r');

%figure (4)
%plot ((1/n).^4,dt);

figure (5)
plot (t,hmin);

figure (6)
plot(t,dt,'d');
xlabel('t [s]');
ylabel('\Deltat [s]');

end
%figure (6)
 %       plot(1./n,hmin)

% figure (7)
% plot(dt.^4,Emec,'bx')
