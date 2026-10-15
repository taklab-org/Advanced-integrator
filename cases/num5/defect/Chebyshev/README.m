clf
% Define a function handle (multi dimension can be accepted)
I = [-1,1]; M = 101;
% f = @(x) [exp(erf(x.^2)+x.^5).*sin(3*pi*x) + x, exp(x)];
f = @(x) exp(erf(x.^2)+x.^5).*sin(3*pi*x) + x;
% f = @(x) cospi(x);
%
% Get the Chebyshev coefficients
% a_{k,l}, |k| < M, l = 1,2,...
a = chebcoeffs(f,M,I);
%
% Plot profile of f(x)
subplot(1,2,1)
plot_cheb(a,I)
%
% Plot
subplot(1,2,2)
plot_chebcoeffs(a)
% hold on, plot_chebcoeffs(a(1:standardChop(a,eps),:))
%
% plot_profile(a)

%% automatic Chebyshev interpolation
clf
% f = @(x) 1./(1+1000*(x+0.5).^2)+1./sqrt(1+1000*(x-0.5).^2);
% f = @(x) exp(erf(x.^2)+x.^5).*sin(3*pi*x) + x;
f = @(x) [exp(erf(x.^2)+x.^5).*sinpi(3*x) + x, exp(x), sinpi(x)];
% f = @(x) sinpi(x);
a = cheb(f,I);
% cheb(@(x) cospi(x))
% cheb(@(x) sinpi(x))
% cheb(@(x) exp(x))
% Plot profile of f(x)
subplot(1,2,1)
plot_cheb(a,I)
%
% Plot
subplot(1,2,2)
plot_chebcoeffs(a)

%% finding roots/endpoints/extrema of Cheyshev interpolant
clf
% f = @(x) 1./(1+1000*(x+0.5).^2)+1./sqrt(1+1000*(x-0.5).^2); I = [-1,1];
% f = @(x) exp(erf(x.^2)+x.^5).*sin(5*pi*x) + x; I = [-1,1];
% f = @(x) exp(x); I = [-1,1];
% f = @(x) besselj(0,x); I = [0,100];
% f = @(x) exp(x).*sech(4*sin(40*x)).^exp(x); I = [-1,1];
% f = @(x)tanh(4*x-1); I = [-1,1];
f = @(x) [exp(erf(x.^2)+x.^5).*sinpi(3*x) + x, exp(x), sinpi(x)];
% f = @(x) [exp(erf(x.^2)+x.^5).*sinpi(3*x) + x, exp(x), sinpi(x), tanh(4*x-1)];
% f = @(x) sinpi(x); I = [-1,1];
%
a = cheb(f,I); ndim = size(a,2);
%
plot_cheb(a,I,length(a)), hold on
%% root of f (red), multi-dimensional available
for i = 1:ndim
    r = chebroots(a(:,i),I);
    plot(r,zeros(length(r)),'o','MarkerSize',8,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6])
end
%
%% endpoints of f (blue), multi-dimensional available
epts = chebendpoints(a);
plot(I(1)*ones(1,length(epts(1,:))), epts(1,:),'o','MarkerSize',8,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
plot(I(2)*ones(1,length(epts(2,:))), epts(2,:),'o','MarkerSize',8,...
    'MarkerEdgeColor','blue',...
    'MarkerFaceColor',[.6 .6 1])
%
% plot(I,exts,'o','MarkerSize',8,...
% 'MarkerEdgeColor','blue',...
% 'MarkerFaceColor',[.6 .6 1])
%% extrema of f (green), multi-dimensional available
for i=1:ndim
    r = chebextrema(a(:,i),I);
    if ~isempty(r)
    fval = f(r);
    plot(r,fval(:,i),'o','MarkerSize',8,...
        'MarkerEdgeColor','green',...
        'MarkerFaceColor',[.6 1 .6])
    end
end
hold off

%% local minima/maxima of Cheyshev interpolant
clf
% f = @(x) exp(erf(x.^2)+x.^5).*sin(5*pi*x) + x; I = [-1,1];
% f = @(x) exp(x); I = [-1,1];
% f = @(x) besselj(0,x); I = [0,100];
f = @(x) [exp(erf(x.^2)+x.^5).*sinpi(3*x) + x, exp(x), sinpi(x)]; I = [-1,1];
% f = @(x) exp(x).*sech(4*sin(40*x)).^exp(x) - 1; I = [0.5,1];
% f = @(x)tanh(4*x-1); I = [-1,1];
a = cheb(f,I);
%
plot_cheb(a,I,length(a)), hold on
%
%%
[fmax,xmax,fmin,xmin] = chebminmax(a,I);
plot(xmax,fmax,'o','MarkerSize',8,...
'MarkerEdgeColor','green',...
'MarkerFaceColor',[.6 1 .6])
%
plot(xmin,fmin,'o','MarkerSize',8,...
'MarkerEdgeColor','red',...
'MarkerFaceColor',[1 .6 .6])
%
hold off

%% Chebyshev interpolant getting from ode functions
clf
% Initial guess by ode45
opts = odeset('AbsTol',1e-12,'RelTol',2.22045e-14);
% % Lotka Volterra
% tspan = [0,8];
% f = ode89(@(t,y) [y(1)-y(1)*y(2); -y(2)+y(1)*y(2)],tspan,[4,1],opts);
% a = odecheb(f,tspan);
% plot_profile(a)
%
% Roesller
tspan = [0,40];
p = [0.2,0.1,2.2];
f = ode89(@(t,y) Rossler(y,p),tspan,[0.1,0.1,0.1],opts);
a = odecheb(f,tspan);
% plot3_profile(a)
%
subplot(1,2,1)
plot_chebcoeffs(a)
subplot(1,2,2)
plot_cheb(a,tspan)
%
%
hold on
% ndim = 1;%size(a,2);
% for i=1:ndim
%     r = chebextrema(a(:,i),tspan);
%     f1 = deval(f,r,i);
%     plot(r,f1,'o','MarkerSize',8,...
%         'MarkerEdgeColor','blue',...
%         'MarkerFaceColor',[.6 .6 1])
% end
% [fmax,xmax,fmin,xmin] = chebminmax(a,tspan);
[fmag,tmax] = chebmag(a,tspan);
for i=1:ndim
    plot(tmax(i),deval(f,tmax(i),i),'o','MarkerSize',8,...
        'MarkerEdgeColor','green',...
        'MarkerFaceColor',[.6 1 .6])
end


%% Chebyshev interpolant getting from dde functions
clf
gamm = 1; beta = 2; tau  = 2;
n = 8.79; %rho
tspan = [0,2];
%
MackeyGlass = @(t,y,z) beta*z(:,1)./(1+z(:,1)^n) - gamm*y;
ddehist = @(t) 0.5;
opts = ddeset('AbsTol',1e-12,'RelTol',2.22045e-14);
% opts = ddeset('RelTol',2.22045e-14);
for i=1:5 % Method of steps
    f = dde23(MackeyGlass,tau,ddehist,tspan,opts);
    %
    a = ddecheb(f,tspan);
    subplot(1,2,1)
    plot_chebcoeffs(a), hold on
    subplot(1,2,2)
    plot_cheb(a,tspan),hold on
    [fmax,xmax] = chebmax(a,tspan);
    plot(xmax,fmax,'o','MarkerSize',8,...
        'MarkerEdgeColor','green',...
        'MarkerFaceColor',[.6 1 .6])
    tspan(1)=tspan(2);
    tspan(2)=tspan(2)+tau;
    ddehist=f;
end

%% Chebyshev interpolant for PDEs in [0,1] getting from ode functions
clf
tspan = [0,2^-4]; N = 15; gamma = exp(1i*(pi/3));
u0 = zeros(2*N-1,1);
u0(N)=50; u0(N-1)=-25; u0(N+1)=-25;
opts = ddeset('AbsTol',1e-12,'RelTol',2.22045e-14);
f = ode89(@(t,y) complex_heat(y,gamma),tspan,u0,opts);
a = pdecheb(f,tspan,3);size(a)
plot_cheb_fourier(a,tspan,[0,1],1)



%% Chebyshev interpolant for higher spatial dimensional PDEs
% TODO



%% ode/pde functions
function dy = Rossler(Y,p)
x = Y(1); y = Y(2); z = Y(3);
a = p(1); b = p(2); c = p(3);
dy = [- y - z;...
    x + a*y;...
    b + x*z - c*z];
end
%
function dy = complex_heat(y,gamma)
N = size(y,1);
fy = quadratic(y,y);
m = (N-1)/2;
k = (-m:m)';
dy = -(4*pi^2)*(k.^2).*y+fy;
dy = dy*gamma;
end
%
function [s]=quadratic(a1,a2)
m=(length(a1)+1)/2;
ta1=[zeros(m,1);a1;zeros(m,1)]; tu1=ifft(ifftshift(ta1));
ta2=[zeros(m,1);a2;zeros(m,1)]; tu2=ifft(ifftshift(ta2));
F=fftshift(fft(tu1.*tu2));
s=(4*m-1)*F(m+1:3*m-1);
end
