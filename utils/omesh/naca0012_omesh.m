clear all
close all
clc

%% define number of grid points and airfoil thickness
%-- nx: number of grid points in circumferential direction: nx must be off
%-- ny: number of grid points in normal direction
%-- thick: NACA symmetric thickness (%)

nx=257;
ny=257
thick=0.12;

if (mod(nx,2)==0)
    error('number of grid points in circumfrential direction must be odd');
end


%% -- Airfoil coordinates based on Cheb points 
NC=nx-1;
m=(nx-1)/2+1;

xTE=1.008930411;
x=linspace(0,xTE,m);

x = xTE/2 + xTE/2 * cos(([1:m]-1)/(m-1)*pi);
y = 5*thick*(0.2969*x.^0.5 -0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);

x=[x x(end-1:-1:1)];
y=[y -y(end-1:-1:1)];
z_naca=x+1i*y;

%%
rLE=1.1019*thick^2;

yxTE=5*thick*(0.2969/2/xTE^0.5 -0.1260 -0.3516*2*xTE +0.2843*3*xTE^2 -0.1015*4*xTE^3);
tau=2*atan(abs(yxTE));
pow=(pi/(2*pi-tau));

z1=xTE;
z2=rLE/2;


eta1=0.77043505;
eta2=0.24642903;
etaLE=0.20139626;
etaTE=eta1;

etaC=0.5*(etaLE+etaTE);

r=((z_naca-z1)./(z_naca-z2)).^pow;
eta_naca = (eta1-r*eta2)./(1-r);
r_naca=abs(eta_naca - etaC);
t_naca=angle(eta_naca-etaC);

%%-- compute radius of concentrtic circles: R(j)
xx=real(eta_naca); a=(max(xx)-min(xx))/2;
yy=imag(eta_naca); b=(max(yy)-min(yy))/2;
dx=xx(2:end)-xx(1:end-1);
dy=yy(2:end)-yy(1:end-1);
s=sum( (dx.^2+dy.^2).^0.5);
R1=s/2/pi;
for j=1:ny
    R(j) = R1*exp((j-1)*(2*pi/NC));
end


%%-- compute R(i,j)
for j=1:ny
    for i=1:nx
        r(i,j)=(r_naca(i)* (R(end)-R1) + R(end)*((R(j)-R1)))/(R(end)-R1);
    end
end


%% compute 2D grid
for j=1:ny;
    for i=1:nx
        eta(i) = etaC + r(i,j)*exp(1i*t_naca(i));
    end
    
    RHS=( (eta-eta1)./(eta-eta2)).^(1/pow);
    
    z=(z1-z2*RHS)./(1-RHS);
    
    X(1:nx,j)=real(z(1:nx));
    Y(1:nx,j)=imag(z(1:nx));
    Y(1 ,j)=(Y(1,j) + Y(nx,j))/2.0;
    Y(nx,j)=Y(1,j);
end 

rpx= abs(X(1,:)-X(nx,:));
rp = ((X(1,:)-X(nx,:)).^2 + (Y(1,:)-Y(nx,:)).^2).^0.5;
max(rpx);

%% enforc symetry
for j=1:ny
    for i=1:(NC/2+1)
        ic=nx+1-i;
        xb(i,j) = (X(i,j)+X(ic,j))/2.0;
        yb(i,j) = (Y(i,j)-Y(ic,j))/2.0;
    end
end

for j=1:ny
    for i=1:(NC/2+1)
        ic=nx+1-i;
        gx(i,j) = xb(i,j);
        gy(i,j) = yb(i,j);
        gx(ic,j)= xb(i,j);
        gy(ic,j)=-yb(i,j);        
    end
end




%%-- plot grid
figure(1);
colormap([1 1 1]);
pc=pcolor(gx(1:nx,1:ny),gy(1:nx,1:ny),0*X(1:nx,1:ny));
pc.LineWidth=1.25;
%pc.Color='b';
xlim([-1.5 1.5]);
ylim([-1.25 1.25]);
axis equal
drawnow
hold on
plot(x,y,'r-',x,-y,'r-','LineWidth',1.35);
hold off;

set(gca,'Visible','off');
fout=['naca0012_omesh'];
% print(gcf,fout,'-dpng','-r800');     

nx1=nx-1;
%% write grid file for fvs2d
ncells=(nx-1)*(ny-1);
nnodes=(nx-1)*ny;
naca=['00' num2str(thick*100)];
fid=fopen(['naca' naca '_omesh.grid'],'w');
fprintf(fid,'%s\n', '#nodes, #triangles_cells, #quad_cells');
fprintf(fid,'%i %i %i\n',nnodes, 0,ncells);
for j=1:ny
    for i=1:nx1
        fprintf(fid,'%e %e\n',gx(i,j),gy(i,j));
    end
end
for j=1:ny-1
    for i=1:nx-2
        ij=(j-1)*nx1+i;
        fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij+nx1+1,ij+1);
    end
    i=nx-1;
        ij =(j-1)*nx1+i;
        ij1=(j-1)*nx1+1;
        fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij1+nx1,ij1);    
end
fclose(fid);


%% write bc file
fid=fopen(['naca' naca '_omesh.bc'],'w');
fprintf(fid,'%s\n','2     !-- #boundaries');
fprintf(fid,'%i %s\n',nx-1,'  slip_wall      !-- boundary 1');
fprintf(fid,'%i %s\n',nx-1,'  freestream     !-- boundary 2');
for j=1:1
    for i=1:nx-1
        ij=(j-1)*nx1+i;
        fprintf(fid,'%i\n',ij);
    end
end
for j=ny-1:ny-1
    for i=1:nx-1
        ij=(j-1)*nx1+i;
        fprintf(fid,'%i\n',ij);
    end
end
fclose(fid);


%% write tecplot file
fid=fopen(['naca' naca '_omesh.plt'],'w');
fprintf(fid,'%s\n', 'variables="x" "y"');
fprintf(fid,'%s %i %s %i %s\n','ZONE NODES=',nnodes, ', ELEMENTS=',ncells, ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL');
for j=1:ny
    for i=1:nx-1
        fprintf(fid,'%e %e\n',gx(i,j),gy(i,j));
    end
end
for j=1:ny-1
    for i=1:nx-2
        ij=(j-1)*nx1+i;
        fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij+nx1+1,ij+1);
    end
    i=nx-1;
    ij =(j-1)*nx1+i;
    ij1=(j-1)*nx1+1;
    fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij1+nx1,ij1);
end
fclose(fid);


%%
for j=1:ny-1
    for i=1:nx-1
        xc(i,j) = (gx(i,j) + gx(i+1,j) + gx(i,j+1) + gx(i+1,j+1))/4.0;
        yc(i,j) = (gy(i,j) + gy(i+1,j) + gy(i,j+1) + gy(i+1,j+1))/4.0;
    end
end

nnodes1=(nx-1)*(ny-1);
ncells1=(nx-1)*(ny-2);

%--
fid=fopen(['naca' naca '_omesh_cellcntr.grid'],'w');
fprintf(fid,'%s\n', '#nodes, #triangles_cells, #quad_cells');
fprintf(fid,'%i %i %i\n',nnodes1, 0,ncells1);
for j=1:ny-1
    for i=1:nx-1
        fprintf(fid,'%e %e\n',xc(i,j),yc(i,j));
    end
end
for j=1:ny-2
    for i=1:nx-2
        ij=(j-1)*nx1+i;
        fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij+nx1+1,ij+1);
    end
    i=nx-1;
    ij =(j-1)*nx1+i;
    ij1=(j-1)*nx1+1;
    fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij1+nx1,ij1);
end
fclose(fid);

%--
fid=fopen(['naca' naca '_omesh_cellcntr.plt'],'w');
fprintf(fid,'%s\n', 'variables="x" "y"');
fprintf(fid,'%s %i %s %i %s\n','ZONE NODES=',nnodes1, ', ELEMENTS=',ncells1, ', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL');
for j=1:ny-1
    for i=1:nx-1
        fprintf(fid,'%e %e\n',xc(i,j),yc(i,j));
    end
end
for j=1:ny-2
    for i=1:nx-2
        ij=(j-1)*nx1+i;
        fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij+nx1+1,ij+1);
    end
    i=nx-1;
    ij =(j-1)*nx1+i;
    ij1=(j-1)*nx1+1;
    fprintf(fid,'%i %i %i %i\n', ij,ij+nx1,ij1+nx1,ij1);
end
fclose(fid);