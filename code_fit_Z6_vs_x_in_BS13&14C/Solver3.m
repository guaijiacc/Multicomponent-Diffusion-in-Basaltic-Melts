
function [X1] = Solver3(D0,a,Cmin,Cmax,x,dx,dt,numx,numx1,numt)
% this subroutine solves a diffusion eq. (for dif couple)
% dC/dt=d(D0*exp(a*C)*dC/dx)/dx; D0 is D at C=0
% using u for C
% Column 1 in X1 (the returned matrix) is x, and column 2 is C 
%    at time=numt*dt s

X1=zeros(numx,3);
ratio=D0*dt/(dx*dx); 
if ratio*exp(a*Cmax)>=0.5
    disp(ratio*exp(a*Cmax)); disp('the algorithm is unstable');
end
u=zeros(numx,1); v=zeros(numx,1); %initialize everything to zero
for i=1:numx1-1
   u(i) = Cmin;
end
for i=numx1+1:numx
    u(i)=Cmax;
end
u(numx1)=0.5*Cmin+0.5*Cmax;
u0=u;

%iterate difference equation
for j=1:numt
    v=(exp(a*u)-1)/a;
    for i=2:numx-1
        u(i)=u(i)+ratio*(v(i+1)+v(i-1)-2*v(i));
    end
    u(1)=(4*u(2)-u(3))/3; u(numx)=(4*u(numx-1)-u(numx-2))/3; 
end

%store x,final conc, and initial conc in matrix X1
X1(:,1)=x; 
X1(:,2)=u;






