function [volume, S, L] = make_laplacian(nx,dx,ny,dy,nz,dz)
%  make_laplacian    calculate 1-D, 2-D or 3-D numerical LaPlacian Operator
% USAGE:  [volume, S, L] = make_laplacian(nx,dx,ny,dy,nz,dz);
% Last modified KCC 4/10/2008 
% calculate the sparse operators of dimension (nx*ny*nz by nx*ny*nz)
% such that 1/volume  integral {Laplacian (model)}^2 dV  = |L*model(:)|^2
% and           1/volume  integral (model)^2 dV  = |S*model(:)|^2
% 
% for 3-D operators, enter all 6 input parameters.
% for 2-D operators, enter only nx,dx,ny,dy or enter nz=1;dz=1;
% for 1-D operator, enter only nx,dx or enter ny=1; nz=1; dy=1;dz=1;
%
% where |v|^2= v'*v
%
% in matrix form the model is organized by model(ix,iy,iz)
% and in vector form as model(:) 
% model is assumed to be parameterized by its value at nodes on a rectangular grid
% with evenly spaced in x, in y and in z
% Input parameters:
%  nx   number of nodes in x direction
%  ny   number of nodes in y direction
%  nz   number of nodes in z direction
%  dx   x-node spacing (km)
%  dy   y-node spacing (km)
%  dz   z-node spacing (km)
%
% Output parameters:
% total volume of grid, and the sparse arrays S and L
%
% In 2-D the actual LaPlacian of Z is:  LZ = reshape(L*Z(:),nx,ny)/sqrt(dx*dy/volume);
% Example 1 below demonstrates this in 2-D for a model of the form Z(x,y)=real( exp(i*(kx*x=ky*y)) ); 
% then the Laplacian of Z is (-kx^2-ky^2)*Z, so if k is a unit vector then the size of Z and the size
% of laplacian Z are the same.  Because  1/T integral 0->T cos(2*pi*t/T)^2 = 0.5, the size of Z 
% by the definitions above and the size of the Laplacian of Z should both be 0.5 (see example below)
% 
% %Example1:
%
% xmax=15; ymax=10; dx=.125; dy=.250; kx=sind(60); ky=cosd(60); x=[0:dx:xmax]';y=[0:dy:ymax]'; nx=length(x); ny=length(y);
% [volume, S, L] = make_laplacian(nx,dx,ny,dy); [X,Y] = ndgrid(x,y); Z=real(exp(sqrt(-1)*(kx*X+ky*Y))); ZL=reshape(L*Z(:),nx,ny); 
% disp(sprintf('model size and roughness are: %f, %f',[norm(S*Z(:)),norm(L*Z(:))].^2))
% figure(2);clf;for kk=1:2; subplot(2,1,kk); if kk==1; surf(X,Y,Z);title('Z');else surf(X,Y,ZL/sqrt(dx*dy/volume)); title('L*Z'); end; 
% view(2);shading('interp');xlabel('x');ylabel('y');axis([-inf inf -inf inf]);set(gca,'dataAspectRatio',[1 1 1]);colorbar; end
%
% %Example2:
%
% nx=5;ny=6;nz=7;dx=1;dy=2;dz=3;  [volume, S, L] = make_laplacian(nx,dx,ny,dy,nz,dz);
% figure(1);clf;spy(L>0,'b');hold on;spy(L<0,'r');title('3-D LaPlacian Operator for 5x6x7 grid')
% figure(2);clf;hist(L(:),200);;set(gca,'ylim',[0 500])
%
% %Example3:
%
% nx=5;ny=6;dx=1;dy=2;  [volume, S, L] = make_laplacian(nx,dx,ny,dy);
% figure(1);clf;spy(L>0,'b');hold on;spy(L<0,'r');title('3-D LaPlacian Operator for 5x6 grid')
% figure(2);clf;hist(L(:),200);;set(gca,'ylim',[0 500])
%
% %Example3:
%
% nx=5;dx=1;ny=4,dy=1;[volume, S, L] = make_laplacian(nx,dx,ny,dy);full(L)*sqrt(length(L))
%
% Features to add:  
%1)  more general options for node spacing (e.g. even spacing in lat
%  and lon, but variable in distance)
% 2) boundary conditions (which could be different at each edge):
%   currently all edges use BC=2 below.
%  BC=0:   model=0; 
%  BC=1:   model'=0; 
%  BC=2:   model''=0; 
%  BC=3;   model=0, model'=0; model''=0 everywhere outside grid
%  BC=4;   minimize curvature for model that wraps around on itself  (e.g. around longitues on earth)
%  each edge could have a different BC
%
% Ken Creager   12/05/00;  modified 4/10/08 kcc@ess.washington.edu

 
BC=2;


if nargin<6; nz=1; dz=1; end
if nargin<4; ny=1; dy=1; end

n    = nx*ny*nz;                % number of nodes
L    = sparse(n,n);             % Laplacian Operator
N    = [nx ny nz]';             % number of x, y, and z nodes
D    = [dx dy dz]';             % x, y, and z node spacing (km)

volume = nx*dx*ny*dy*nz*dz;  % total volume of all nodes

 dV   = sqrt( dx*dy*dz/ volume);  % sqrt of dV/V for each node ( same as 1/sqrt(n) ) 

 ndim=nargin/2;   % number of dimensions
 
 if     ndim==1; [ix]       = [1:nx]';  
 elseif ndim==2; [ix,iy]    = ndgrid( 1:nx, 1:ny );  
 elseif ndim==3; [ix,iy,iz] = ndgrid( 1:nx, 1:ny, 1:nz );  
 end
 
 L     = sparse(n,n);             % Laplacian Operator

 for k=1:ndim;
   if k==1;     I=ix; dd=dx; nn=nx;
   elseif k==2; I=iy; dd=dy; nn=ny;
   elseif k==3; I=iz; dd=dz; nn=nz;
   end
   i0=find(I(:)>1 & I(:)<nn); 
   ip=find(I(:)>2);
   im=find(I(:)<nn-1);
   dV_d1sq = dV/(dd^2);
   L(i0+(im-1)*n) =                    dV_d1sq;
   L(i0+(ip-1)*n) =                    dV_d1sq;
   L(i0+(i0-1)*n) = L(i0+(i0-1)*n) - 2*dV_d1sq;
 end
 
S       = spdiags(sqrt(1/n)+zeros(n,1),0,n,n);

return
 



% loop over x, skipping first and last rows
% boundary condition is to assume Laplacian along all edge nodes = 0

%x = primary, y/z secondary (I=123);  x=1; y=2; z=3;

%I  = [primary dimension; 2 secondary dimensions in order] = xyz; yxz; zxy or 123; 213 312
%J  = [123  213  231


%1=x; 2=y; 3=z; I=123(xyz);  x=1; y=2; z=3; J=123;
%1=y; 2=x; 3=z; I=213(yxz);  x=2; y=1; z=3; J=213
%1=z; 2=x; 3=y; I=312(zxy);  x=2; y=3; z=1; J=231;

%1=z; 2=y; 3=x; I=321;  x=3; y=2; z=1; J=321;
for ii=1:3                  % loop over the 3 dimensions
 if     ii==1; I = [1 2 3]; % find 2nd-deriviteve in x coordinate; looping over [yz] plane; index for xyz is 123
 elseif ii==2; I = [2 1 3]; % find 2nd-deriviteve in y coordinate; looping over [xz] plane; index for yxz is 213
 elseif ii==3; I = [3 1 2]; % find 2nd-deriviteve in z coordinate; looping over [xy] plane; index for zxy is 312
 end

 % primary coordinate is the one along which the derivative is being taken
 
 d1          =D(I(1)); % node spacing along primary coordinate 
 d2          =D(I(2)); % node spacing  along one secondary coordinate
 d3          =D(I(3)); % node spacing  along other secondary coordinate
 
 % get index vectors into place corresponding to secondary corrdinates ;  must have I(2) < I(3) to keep order in matrix
 [ind2,ind3] = ndgrid( 1:N(I(2)) , 1:N(I(3)) );  
 n1          = N(I(1));                      % number of points in primary coordinate
 n23         = N(I(2)) * N(I(3));            % number of points in secondary coordinate plane
 ind         = zeros(n23,3);                 % index into primary coordinate
 ind(:,2)    = ind2(:);                      % index into one secondary corrdinate
 ind(:,3)    = ind3(:);                      % index into other secondary coordinate
 
 indx = ind(:,find(I==1));                   % index into x-coordinate
 indy = ind(:,find(I==2));                   % index into y-coordinate
 indz = ind(:,find(I==3));                   % index into z-coordinate
 
 ix = (ii==1);                               % 1 if x is primary coordinate , otherwize =0
 iy = (ii==2);                               % 1 if y is primary coordinate , otherwize =0
 iz = (ii==3);                               % 1 if z is primary coordinate , otherwize =0
 if n1>1;
 for j = 1:n1;                               % loop over primary coordinate (avoiding edges, boundary condition is to minimize second derivative everywhere)
  j0=j+0; jm=j-1; jp=j+1;                     % index into primary coordinate for point of interest and  the two adjacent points
  i0=sub2ind([nx,ny,nz],indx+ix*j0,indy+iy*j0,indz+iz*j0); % get index into model vector for points in this plane
  if j==1; im=[];else
  im=sub2ind([nx,ny,nz],indx+ix*jm,indy+iy*jm,indz+iz*jm); % get index into model vector for points in adjacent plane
  end 
  if j==n1; ip=[];else
  ip=sub2ind([nx,ny,nz],indx+ix*jp,indy+iy*jp,indz+iz*jp); % get index into model vector for points in adjacent plane
  end
  dV   = sqrt( d1*d2*d3/ volume);                          % scaled volume of each node
  % square root of fractional volume for these nodes
  % square root is because integral Laplacian^2 dV  =  sum over nodes [Laplacian  sqrt(dA)]^2;  dV is dimensionless
  % divide by total volume so that the result is normalized by volume and has units of roughness is  (model / length^2) ^2
  dV_d1sq = dV/(d1^2);
  if j>1 & j<n1;
   L(i0+(im-1)*n) =                    dV_d1sq;
   L(i0+(ip-1)*n) =                    dV_d1sq;
   L(i0+(i0-1)*n) =L(i0+(i0-1)*n) - 2* dV_d1sq;
  
  % for k=1:n23;                          % loop over all points in secondary plane
  %   J=i0(k);
  %   L(J,ip(k)) =            dV_d1sq;
  %   L(J,im(k)) =            dV_d1sq;
  %   L(J,J)     = L(J,J) - 2*dV_d1sq;
  % end
  elseif BC==2;
  elseif BC==3;% keyboard
   if j==1;
    for k=1:n23;                          % loop over all points in secondary plane
      J=i0(k);
      L(J,ip(k)) =            dV_d1sq;
      L(J,J)     = L(J,J) - 2*dV_d1sq;
    end
   elseif j==n1;
    for k=1:n23;                          % loop over all points in secondary plane
      J=i0(k);
      L(J,im(k)) =            dV_d1sq;
      L(J,J)     = L(J,J) - 2*dV_d1sq;
    end
   end
  end
 end
 end
end

% now construct operator for making the model small
ii      = [1:n]';
[i,j,k] = ind2sub([nx,ny,nz],ii);
% calculate volume of each cell using the indices i,j,k
dV      = sqrt(dx*dy*dz/volume);  dV=dV+zeros(size(ii));
S       = spdiags(dV,0,n,n);

return


% OLD CODE FOR 2-D problem 

[yind,zind] = ndgrid(1:ny,1:nz); yind=yind(:);zind=zind(:);xind=zeros(size(yind));
for kx=2:nx-1;
  i0=sub2ind([nx,ny,nz],(kx+0)+xind,yind,zind); % get index into model vector for points in this slice
  im=sub2ind([nx,ny,nz],(kx-1)+xind,yind,zind); % get index into model vector for  adjacent points to the south
  ip=sub2ind([nx,ny,nz],(kx+1)+xind,yind,zind); % get index into model vector for  adjacent points to the north
  dV   = sqrt(dx*dy*dz / volume);
  % square root of fractional volume for these nodes
  % square root is because integral Laplacian^2 dV  =  sum over nodes [Laplacian  sqrt(dA)]^2;  dV is dimensionless
  % divide by total volume so that the result is normalized by volume and has units of roughness is  (model / length^2) ^2
  tmp  = dV./(dx.^2); 
  if length(tmp)==1; tmp = tmp+zeros(size(xind)); end
  for k=1:length(xind);                     % loop over all longitudes
    I=i0(k);
    L(I,ip(k)) =            tmp(k);
    L(I,im(k)) =            tmp(k);
    L(I,I)     = L(I,I) - 2*tmp(k);
  end
end
 
oney = ones(ny,1);
yvec = [1:ny]';
for kx=2:nx-1;         % loop over each column (longitude)  skipping first and last columns
  i0=sub2ind([ny,nx],yvec,(kx+0)*oney); % get index into model vector for points in this row
  im=sub2ind([ny,nx],yvec,(kx-1)*oney); % get index into model vector for adjacent points to the south
  ip=sub2ind([ny,nx],yvec,(kx+1)*oney); % get index into model vector for adjacent points to the north
  dA   = sqrt(dy(kx)*dx / volume);
  tmp  = dA./(dx.^2);
  for k=1:ny;          % loop over all latitudes
    I=i0(k);
    L(I,ip(k)) =            tmp(k);
    L(I,im(k)) =            tmp(k);
    L(I,I)     = L(I,I) - 2*tmp(k);
  end
end

% now construct operator for making the model small
onex = ones(nx,1);
xvec = [1:nx]';
Sdiag= zeros(n,1);
for ky=1:ny;
  i    = sub2ind([ny,nx],(ky+0)*onex,xvec); % get index into model vector for points in this row
  dA   = sqrt(dx(ky)*dy / volume);
  Sdiag(i,1) = dA(:);
end
S      = spdiags(Sdiag,0,n,n);


return

%%%%%%%%%%%%%%%%%%%%
%
% OLD CODE  11/16/2000
%
%%%%%%%%%%%%%%%%%%%%

for kk=1:2;
  if kk==1; xx=lon_in;  yy=lat_in;  XX=LON_IN;  YY=LAT_IN; 
  else;     xx=lon_out; yy=lat_out; XX=LON_OUT; YY=LAT_OUT;
  end
  nyy=length(yy);     % number of nodes in y direction
  nxx=length(xx);     % number of nodes in x direction
  nn =nxx*nyy;        % number of nodes
  RR =sparse(nn,nn);
  SS =zeros(nn,1);
  dy = R*y_step*pi/180 * ones(size(xx));  % y spacing (km) for each x-node
  dx = R*x_step*pi/180 * cos(yy*pi/180);  % x spacing (km) for each y-node
%  dx = R*x_step*pi/180 * ones(size(yy));  % x spacing (km) for each y-node TEST CODE
  area = sum(sum(dx(:)*dy(:)'));          % total area of all nodes (km^2)

  % loop over each row (latitude), skipping first and last rows 
  % boundary condition is to assume lapalcian along all edge nodes = 0
  for ky=2:nyy-1; 
    y=yy(ky);             % get y value for this row
    i=find(YY==y);        % get index into model vector for points in this row
    j=find(YY==y+y_step); % get index into model vector for adjacent points to the north
    k=find(YY==y-y_step); % get index into model vector for adjacent points to the south
    dA   = sqrt(dx(ky)*dy / area); % square root of fractional area for these nodes
    % square root is because integral Laplacian^2 dA  =  sum over nodes [Laplacian  sqrt(dA)]^2;  dA is dimensionless
    % divide by total area so that the result is normalized by area and has units of roughness is  (model / length^2) ^2
    tmp  = dA./(dy.^2);            
    for l=1:nxx;                     % loop over all longitudes
      il=i(l); jl=j(l); kl=k(l);
      RR(il,jl)=tmp(l);
      RR(il,kl)=tmp(l);
      RR(il,il)=RR(il,il)-2*tmp(l);
    end
  end
  
  for kx=2:nxx-1;         % loop over each column (longitude)  skipping first and last columns
    x=xx(kx);             % get x value for this column
    i=find(XX==x);        % get index into model vector for points in this column
    j=find(XX==x+x_step); % get index into model vetor for adjacent points to the east
    k=find(XX==x-x_step); % get index into model vetor for adjacent points to the west
    dA   = sqrt(dy(kx)*dx / area);
    tmp  = dA./(dx.^2);
    for l=1:nyy;          % loop over all latitudes
      il=i(l); jl=j(l); kl=k(l);
      RR(il,jl)=tmp(l);
      RR(il,kl)=tmp(l);
      RR(il,il)=RR(il,il)-2*tmp(l);
    end
  end
  if kk==1; Rin0 =RR;  areain =area;
  else;     Rout0=RR;  areaout=area;
  end
  
% now make the model small  
  SSdiag = zeros(nn,1);
  for ky=1:nyy;
    y=yy(ky);      % get y value for this row
	i=find(YY==y);        % get index into model vector for points in this row
    dA   = sqrt(dx(ky)*dy / area);
	SSdiag(i,1) = dA(:);
  end
  if kk==1; Sin0  = spdiags(SSdiag,0,nn,nn);
  else;     Sout0 = spdiags(SSdiag,0,nn,nn);
  end
end
clear RR
