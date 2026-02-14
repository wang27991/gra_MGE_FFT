clc
clear
close all
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nz=50;
II0=60;
DD0=10;
II=60;
DD=10;
%%
zmin=-3000;
%%
bbb='30_20_10_inv';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xyzg=load(['grav' bbb '.txt'],'-ascii');
xyh0=load(['topo' bbb '.txt'],'-ascii');
outg=['GgF' bbb '.mat'];
nobs=length(xyzg(:,1))-1;
xmin=min(xyh0(2:end,1));
xmax=max(xyh0(2:end,1));
ymin=min(xyh0(2:end,2));
ymax=max(xyh0(2:end,2));
zmax=max(xyh0(2:end,3));
nx=xyh0(1,1);
ny=xyh0(1,2);
N=nx*ny*nz;
%%
lx=(xmax-xmin)/(nx-1);
ly=(ymax-ymin)/(ny-1);
lz=(zmax-zmin)/(nz);
size0=[lx,ly,lz];
x=(xmin-lx/2):lx:(xmax+lx/2);
y=(ymin-ly/2):ly:(ymax+ly/2);
z=(zmax):-lz:(zmin);
GgT0=zeros(nx+1,ny+1,nz+1);
GgT=zeros(2*nx-1,2*ny-1,nz);
tic

xyz0=xyzg(2,1:3);


if nx == ny
    
    for k=1:nz+1
        
        xyz00 = xyz0 - [x(1), y(1), z(k)];
        x00 = xyz00(1); y00 = xyz00(2); z00 = xyz00(3);
        r = sqrt(sum(xyz00.^2));
        GgT0(1,1,k) = y00*log(x00+r) + x00*log(y00+r) - z00*atan(x00*y00/(z00*r));

        xyz00 = xyz0 - [x(1), y(2), z(k)];
        x00 = xyz00(1); y00 = xyz00(2); z00 = xyz00(3);
        r = sqrt(sum(xyz00.^2));
        GgT0(1,2,k) = y00*log(x00+r) + x00*log(y00+r) - z00*atan(x00*y00/(z00*r));

        GgT0(2,1,k) = GgT0(1,2,k);
        GgT0(2,2,k) = -GgT0(1,1,k) - GgT0(1,2,k) - GgT0(2,1,k);
        
        for i=1:nx+1
            for j=i:ny+1  

                if (i==1 && j==1) || (i==1 && j==2) || (i==2 && j==2)
                    continue;
                end
                
                xyz00=xyz0-[x(i),y(j),z(k)];
                x00=xyz00(1); y00=xyz00(2); z00=xyz00(3);
                r=sqrt(sum(xyz00.^2));
                
                value = y00*log(x00+r) + x00*log(y00+r) - z00*atan(x00*y00/(z00*r));
                GgT0(i,j,k) = value;
                
                if i ~= j
                    GgT0(j,i,k) = value;
                end
            end
        end
    end
else

    for k=1:nz+1

        
        min_n = min(nx+1, ny+1);
        

        for i=1:min_n
            for j=i:min_n  
                xyz00 = xyz0 - [x(i), y(j), z(k)];
                x00 = xyz00(1); y00 = xyz00(2); z00 = xyz00(3);
                r = sqrt(sum(xyz00.^2));
                
                value = y00*log(x00+r) + x00*log(y00+r) - z00*atan(x00*y00/(z00*r));
                GgT0(i,j,k) = value;
                
                if i ~= j
                    GgT0(j,i,k) = value;
                end
            end
        end
        

        if nx+1 > min_n  
            for i=min_n+1:nx+1
                for j=1:ny+1
                    xyz00 = xyz0 - [x(i), y(j), z(k)];
                    x00 = xyz00(1); y00 = xyz00(2); z00 = xyz00(3);
                    r = sqrt(sum(xyz00.^2));
                    value = y00*log(x00+r) + x00*log(y00+r) - z00*atan(x00*y00/(z00*r));
                    GgT0(i,j,k) = value;
                end
            end
        elseif ny+1 > min_n  
            for j=min_n+1:ny+1
                for i=1:nx+1
                    xyz00 = xyz0 - [x(i), y(j), z(k)];
                    x00 = xyz00(1); y00 = xyz00(2); z00 = xyz00(3);
                    r = sqrt(sum(xyz00.^2));
                    value = y00*log(x00+r) + x00*log(y00+r) - z00*atan(x00*y00/(z00*r));
                    GgT0(i,j,k) = value;
                end
            end
        end
    end
end

GgT0 = -1000 * (6.67*10^-6) * GgT0; 

GgT00 = -GgT0(1:nx, 1:ny, 1:nz) - GgT0(2:nx+1, 2:ny+1, 1:nz)...
    - GgT0(2:nx+1, 1:ny, 2:nz+1) - GgT0(1:nx, 2:ny+1, 2:nz+1)...
    + GgT0(2:nx+1, 1:ny, 1:nz) + GgT0(1:nx, 2:ny+1, 1:nz)...
    + GgT0(1:nx, 1:ny, 2:nz+1) + GgT0(2:nx+1, 2:ny+1, 2:nz+1);


GgT(nx:2*nx-1, ny:2*ny-1, :) = GgT00;


if ny > 1
    GgT(nx:2*nx-1, 1:ny-1, :) = fliplr(GgT00(:, 2:end, :));
end


if nx > 1
    GgT(1:nx-1, ny:2*ny-1, :) = flipud(GgT00(2:end, :, :));
end


if nx > 1 && ny > 1
    GgT(1:nx-1, 1:ny-1, :) = rot90(GgT00(2:end, 2:end, :), 2);
end

GgF = fft2(GgT);
toc

save(outg, 'GgF');