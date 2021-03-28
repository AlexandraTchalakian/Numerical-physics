

G=6.674e-11;
mt=5,9736e+24;
v1=5000;
r0=600000000;
alpha=1.5*pi/180.0;
ml=7.2377e+22;
RT=6378.1*10e+3;


r1=(-2*G*mt-sqrt(4*(G*mt)^2+4*(v1^2-2*G*mt/r0)*(r0*v1*sin(alpha))^2))/(2*(v1^2-2*G*mt/r0));
hmin1=RT+r1
r2=(-2*G*mt+sqrt(4*(G*mt)^2+4*(v1^2-2*G*mt/r0)*(r0*v1*sin(alpha))^2))/(2*(v1^2-2*G*mt/r0));
hmin2=RT+r2
%  C=[omega^2,
%     -2*omega^2*(x3+x4),
%     omega^2*(x3^2+x2^2+4*x2x3),
%     G*(mt+ml)-2*omega^2*(x2*x3^2+x3*x2^2),
%     omega^2*x2^2*x3^2-2*G*(x3*mt+x2*ml),
%     G*(mt*x3^2+ml*x2^2)
%     ]