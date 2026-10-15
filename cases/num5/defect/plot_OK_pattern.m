function plot_DC_pattern(a,L) % input is one-sided (cosine) Fourier
[m,n] = size(a);
L1 = L(1); L2 = L(2);
%
m_pad = 200;
n_pad = 200;
%
Nx = m + m_pad;
Ny = n + n_pad; % size of consine
%
a_pad = zeros(Nx,Ny);
a_pad(1:m,1:n) = a;
%
dx = pi/(Nx-1)/L1;
dy = pi/(Ny-1)/L2;
%
x = dx*(0:Nx-1);
y = dy*(0:Ny-1);
%
%
aa = zeros(2*Nx-2,2*Ny-2);
%
Sx = Nx;
Sy = Ny;
%
aa(Sx:Sx+Nx-2,Sy:Sy+Ny-2)=a_pad(1:end-1,1:end-1);
aa(Sx-Nx+1:Sx,Sy:Sy+Ny-2)=flip(a_pad(:,1:end-1),1);
aa(Sx:Sx+Nx-2,Sy-Ny+1:Sy)=flip(a_pad(1:end-1,:),2);
aa(Sx-Nx+1:Sx,Sy-Ny+1:Sy)=flip(flip(a_pad,1),2);
u = (2*Nx-2)*(2*Ny-2)*real(ifft2(ifftshift(aa)));
% mesh(x,y,u(Sx:end,Sy:end))
surf(x,y,permute(u(1:Sx,1:Sy),[2,1]),'EdgeColor','none'), hold on
surf(-x,y,permute(u(1:Sx,1:Sy),[2,1]),'EdgeColor','none')
surf(x,-y,permute(u(1:Sx,1:Sy),[2,1]),'EdgeColor','none')
surf(-x,-y,permute(u(1:Sx,1:Sy),[2,1]),'EdgeColor','none')
%
surf(x+pi/L1,y,permute(u(Sx:-1:1,1:Sy),[2,1]),'EdgeColor','none')
surf(x+pi/L1,-y,permute(u(Sx:-1:1,1:Sy),[2,1]),'EdgeColor','none')
surf(-x-pi/L1,y,permute(u(Sx:-1:1,1:Sy),[2,1]),'EdgeColor','none')
surf(-x-pi/L1,-y,permute(u(Sx:-1:1,1:Sy),[2,1]),'EdgeColor','none')
%
surf(x,y+pi/L2,permute(u(1:Sx,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(-x,y+pi/L2,permute(u(1:Sx,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(x,-y-pi/L2,permute(u(1:Sx,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(-x,-y-pi/L2,permute(u(1:Sx,Sy:-1:1),[2,1]),'EdgeColor','none')
%
surf(x+pi/L1,y+pi/L2,permute(u(Sx:-1:1,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(x+pi/L1,-y-pi/L2,permute(u(Sx:-1:1,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(-x-pi/L1,y+pi/L2,permute(u(Sx:-1:1,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(-x-pi/L1,-y-pi/L2,permute(u(Sx:-1:1,Sy:-1:1),[2,1]),'EdgeColor','none')
%
surf(x,y+2*pi/L2,permute(u(1:Sx,1:Sy),[2,1]),'EdgeColor','none')
surf(-x,y+2*pi/L2,permute(u(1:Sx,1:Sy),[2,1]),'EdgeColor','none')
surf(x+pi/L1,y+2*pi/L2,permute(u(Sx:-1:1,1:Sy),[2,1]),'EdgeColor','none')
surf(-x-pi/L1,y+2*pi/L2,permute(u(Sx:-1:1,1:Sy),[2,1]),'EdgeColor','none')
%
surf(x,y+3*pi/L2,permute(u(1:Sx,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(-x,y+3*pi/L2,permute(u(1:Sx,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(x+pi/L1,y+3*pi/L2,permute(u(Sx:-1:1,Sy:-1:1),[2,1]),'EdgeColor','none')
surf(-x-pi/L1,y+3*pi/L2,permute(u(Sx:-1:1,Sy:-1:1),[2,1]),'EdgeColor','none'),hold off
% hold on
% xlabel('$x_1$','interpreter','latex', 'FontSize', 30)
% ylabel('$x_2$','interpreter', 'latex', 'FontSize', 30)
% zlabel('$u$','interpreter', 'latex', 'FontSize', 30)
xticks([]), yticks([])
set(gca,'FontSize',20)
axis tight
