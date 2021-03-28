format LONGG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define some variables
workingfolder = './'; % Path to the folder that contains the binary 
% of the code and the simulations (must end with /). Here we execute the
% script in the same folder.
binfilename = 'Gravitation'; % Name of the binary executable
ndigit = 18; % Precision used in numerical to string conversion for input data file name
% Define parameters: scanned and constant ones
% NB: scanned parameter must always be the first in these lines.

%here, scan on dt, other parameters are kept constant.
%lambda = linspace(0,90,1);
lambda=180;
delta=[0,300e+3,1000e+3];
%epsilon = linspace(100,100,1);
nsimul = length(delta); %number of points in the scan.

dt = repmat(500,1,nsimul);
%tfin = repmat(2419200,1,nsimul); %28jours
tfin = repmat(31536000,1,nsimul); %1 an
G = repmat(6.674e-11,1,nsimul);
Mt = repmat(5.9736e+24,1,nsimul);
Ma = repmat(1e+5,1,nsimul);
Ml = repmat(7.3477e+22,1,nsimul);
Rt = repmat(6378.1e+3,1,nsimul);
Rl = repmat(1738e+3,1,nsimul);
alpha = repmat(90*(pi/180),1,nsimul);
epsilon = repmat(100,1,nsimul);

d = 384748e+3;
rt=-Ml(1)*d/(Mt(1)+Ml(1));
rl=Mt(1)*d/(Mt(1)+Ml(1));
omega=G(1)*(Ml(1)+Mt(1))/d^3;

% C=[omega,
%    -2*omega*(rl+rt),
%    omega*(rl^2+rt^2+4*rt*rl),
%    G*(Mt(1)+Ml(1))-2*omega*(rt*rl^2+rl*rt^2),
%    omega*rl^2*rt^2-2*G(1)*(rl*Mt(1)+rt*Ml(1)),
%    G(1)*(Mt(1)*rl^2+Ml(1))*rt^2];
% z=roots(C);
L3=-386695881.90398;

% ra0x = 6e+8.*cos(lambda*(pi/180));
% ra0y = 6e+8.*sin(lambda*(pi/180));
ra0x = repmat(L3,1,nsimul);
ra0y = delta;
%rt0x = repmat(0,1,nsimul);
%rt0y = repmat(0,1,nsimul);
%rl0x = repmat(384748e+3,1,nsimul);
%rl0y= repmat(0,1,nsimul);

%va0x = repmat(-5e+3*cos(1.4*((2*pi)/360)),1,nsimul);
%va0y = repmat(5e+3*sin(1.4*((2*pi)/360)),1,nsimul);
va0x = repmat(0,1,nsimul);
va0y = -sqrt(omega)*sqrt(ra0x.^2+ra0y.^2)

%va0y = repmat(0,1,nsimul);
%vt0x = repmat(0,1,nsimul);
%vt0y = repmat(0,1,nsimul);
%vl0x = repmat(0,1,nsimul);
%vl0y= repmat(0,1,nsimul);

%alpha_c = (sqrt(((5e+3*5e+3)-2*G(1)*(Mt)/(ra0x(1)))*Rt*Rt+2*G(1)*Mt(1)*Rt(1)))/(ra0x(1)*ra0x(1)*5e+3*5e+3);

% va0x = -5e+3*(cos(lambda*(pi/180)- alpha));
% va0y = -5e+3*(sin(lambda*(pi/180)- alpha)); 
adapt = repmat(true,1,nsimul);

%Exercice 3
rt0x = repmat(-(Ml(1)*d)/(Ml(1)+Mt(1)),1,nsimul);
rt0y = repmat(0,1,nsimul);
rl0x = repmat((Mt(1)*d)/(Ml(1)+Mt(1)),1,nsimul); 
rl0y = repmat(0,1,nsimul);

%vt0x = repmat(-Ml(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*sin(1.4*((2*pi)/360)),1,nsimul);
%vt0y = repmat(Ml(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*cos(1.4*((2*pi)/360)),1,nsimul);
%vl0x = repmat(Mt(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*sin(1.4*((2*pi)/360)),1,nsimul);
%vl0y = repmat(-Mt(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1))))*cos(1.4*((2*pi)/360)),1,nsimul);

vt0x = repmat(0,1,nsimul);
vt0y = repmat(-Ml(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1)))),1,nsimul);
vl0x = repmat(0,1,nsimul);
vl0y = repmat(Mt(1)*d*sqrt((G(1))/(d^3*(Mt(1)+Ml(1)))),1,nsimul);


%vt0x = repmat ();
save('w.mat')
infilename = { 'delta' }; % Variable name whose value will be included in the simulation file name
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
n = [];
hmin = [];
Xa_2 =[];
Ya_2 = [];



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
h = data(:,15);
Emec = data(:,16);
hla = data(:,17);

% Xaxx = zeros(length(Xa));
% Yayy = zeros(length(Ya));
% Xlxx = zeros(length(Ya));
% Ylyy = zeros(length(Ya));
% Xtxx = zeros(length(Ya));
% Ytyy = zeros(length(Ya));

% for i = 1 : length(Xa)
% Xaxx(i) = Xa(i)*cos(sqrt(omega)*t(i))+Ya(i)*sin(sqrt(omega)*t(i));
% Yayy(i)=-Xa(i)*sin(sqrt(omega)*t(i))+Ya(i)*cos(sqrt(omega)*t(i));
% Xlxx(i)=Xl(i)*cos(sqrt(omega)*t(i))+Yl(i)*sin(sqrt(omega)*t(i));
% Ylyy(i)=-Xl(i)*sin(sqrt(omega)*(i))+Yl(i)*cos(sqrt(omega)*t(i));
% Xtxx(i)=Xt(i)*cos(sqrt(omega)*t(i))+Yt(i)*sin(sqrt(omega)*t(i));
% Ytyy(i)=-Xt(i)*sin(sqrt(omega)*t(i))+Yt(i)*cos(sqrt(omega)*t(i));
% end

Xaxx=Xa.*cos(sqrt(omega)*t)+Ya.*sin(sqrt(omega)*t);
Yayy=-Xa.*sin(sqrt(omega)*t)+Ya.*cos(sqrt(omega)*t);
Xlxx=Xl.*cos(sqrt(omega)*t)+Yl.*sin(sqrt(omega)*t);
Ylyy=-Xl.*sin(sqrt(omega)*t)+Yl.*cos(sqrt(omega)*t);
Xtxx=Xt.*cos(sqrt(omega)*t)+Yt.*sin(sqrt(omega)*t);
Ytyy=-Xt.*sin(sqrt(omega)*t)+Yt.*cos(sqrt(omega)*t);

jj = 1;

% while and(and(h(jj) - Rt(1) > 0,  jj <= length(t)),hla(jj) - Rl(1) > 0)
% 
%     Xa_2(jj) = Xa(jj);
%     Ya_2(jj) = Ya(jj);
%     jj = jj+1;
%     jj;
% end



t2 = linspace(1,10,1000);

n(ii) = length(dt);

figure (1)
hold on;
p1 = plot(Xlxx,Ylyy,'-');
hold on;
p2 = plot(Xtxx,Ytyy,'-');
hold on;
%p3 = plot(Rl(1)*cos(t2)+Xl(1),Rl(1)*sin(t2)+Yl(1));
hold on;
%p4 = plot(Rt(1)*cos(t2)+Xt(end),Rt(1)*sin(t2)+Yt(end));
hold on;
p5=plot(Xaxx,Yayy,'-');
xlabel('X [m]');
ylabel('Y [m]');
legend([p1(1),p2(1),p5(1)],'trajectoire lune','trajectoire Terre','trajectoire astéroïde');
axis equal;

figure (2)
hold on;
t1 = plot(Xl,Yl,'-');
hold on;
t2 = plot(Xt,Yt,'-');
hold on;
t3=plot(Xa,Ya,'--');
axis equal
xlabel('X [m]')
ylabel('Y [m]')
legend([t1(1),t2(1),t3(1)],'trajectoire lune','trajectoire terre','trajectoire astéroïde')
figure (3)
hold on;
plot(t,dt,'d')

figure (5)
hold on;
plot(Xa_2,Ya_2);
hold on;
p1 = plot(Xl,Yl,'g');
hold on;
p2 = plot(Xt,Yt,'r');

end
