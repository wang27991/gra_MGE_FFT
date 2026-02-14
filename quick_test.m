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


%%
zmin=-3000;
Itermax=10;
%%
bbb='30_20_10_inv';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['GgF' bbb '.mat']);
% EMEKY=load("777.dat");
% AMEKY=repmat(100,1,2400)';
% QMEKY=EMEKY(1:2400,1);
% CMEKY=EMEKY(1:2400,2);
% WMEKY=EMEKY(1:2400,3);
% XMEKY=[QMEKY,CMEKY,AMEKY,WMEKY];
% kkr=[40,60,0,0];
% xyzg=[kkr;XMEKY];
xyzg=load(['grav' bbb '.txt'],'-ascii');
xyh0=load(['topo' bbb '.txt'],'-ascii');
outginv=['gravinv' bbb '.txt'];
%%
g=xyzg(2:end,4);
xmin=min(xyzg(2:end,1));
xmax=max(xyzg(2:end,1));
ymin=min(xyzg(2:end,2));
ymax=max(xyzg(2:end,2));
zmax=max(xyh0(2:end,3));
nx=xyh0(1,1);%%网格
ny=xyh0(1,2);
N=nx*ny*nz;
lx=(xmax-xmin)/(nx-1);
ly=(ymax-ymin)/(ny-1);
lz=(zmax-zmin)/(nz);
x=(xmin):lx:(xmax);
y=(ymin):ly:(ymax);
z=(zmax-lz/2):-lz:(zmin+lz/2);
[xxx,yyy,zzz]=meshgrid(x,y,z);
%%
nx0=xyzg(1,1); ny0=xyzg(1,2);
tp=reshape(g,ny0,nx0);
Wzf0=(zmax-z+abs(zmin)*0.01);
%%
Wzfg=Wzf0.^1.2;
%%
tic
figure(1)
Cg=inversionGrF(GgF,g,Itermax,Wzfg);%稀疏的Gr
toc
%% 拟合差
% rho=Cg;
% r0=Gg*Cg-g;
% display(sqrt(r0'*r0/length(g)));
% m=Cm;
% r0=Gm*Cm-T;
% display(sqrt(r0'*r0/length(T)));
%%
% figure(2)
% nx0=xyzg(1,1); ny0=xyzg(1,2);
% gg1=reshape(g,ny0,nx0);
% subplot(321)
% imagesc(gg1)
% colorbar
% colormap jet
% TT1=reshape(T,ny0,nx0);
% subplot(322)
% imagesc(TT1)
% colorbar
% colormap jet
% gg2=reshape(Gg*Cg,ny0,nx0);
% subplot(323)
% imagesc(gg2)
% colorbar
% colormap jet
% TT2=reshape(Gm*Cm,ny0,nx0);
% subplot(324)
% imagesc(TT2)
% colorbar
% colormap jet
% subplot(325)
% imagesc(gg1-gg2)
% colorbar
% colormap jet
% subplot(326)
% imagesc(TT1-TT2)
% colorbar
% colormap jet
%% g
Rat=[xmax-xmin,ymax-ymin,zmax-zmin];%设置图片显示比例
rrr=permute(Cg,[2 1 3]);

figure(3)

d=slice(xxx,yyy,zzz,rrr,[],y(round(ny/2)),[]);
% %d=slice(xxx,yyy,zzz,rrr,[],y(round(ny/2)),[]);
% d=slice(xxx,yyy,zzz,rrr,[],[],-1000);


set(d,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
load('mod1.mat');
ck=[0,1,1];
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
set(gca,'PlotBoxAspectRatio',Rat);
%% 保存
xyzdg=[xxx(:),yyy(:),zzz(:),rrr(:)];
save(outginv,'xyzdg','-ascii');
