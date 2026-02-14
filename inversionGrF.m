function Cg1=inversionGrF(Gg,g,Imax,Wzg1)
num=0;
k=0;
[nx,ny,nz]=size(Gg);
nx=(nx+1)/2;
ny=(ny+1)/2;

Cg0=zeros(2*nx-1,2*ny-1,nz);
r1g0=zeros(2*nx-1,2*ny-1,nz);
p1g=zeros(2*nx-1,2*ny-1,nz);
Cg1=zeros(nx,ny,nz);

g0=reshape(g,ny,nx);
g1=g0';

% g1=zeros(2*nx-1,2*ny-1); T1=zeros(2*nx-1,2*ny-1);
g00=zeros(2*nx-1,2*ny-1);
% g1(1:nx,1:ny)=g0;        T1(1:nx,1:ny)=T0;

Vg=Gg;

for i=1:nz
    Vg(:,:,i)=Gg(:,:,i)*Wzg1(i);
end

rrg=zeros(Imax,1);

ag=1e-1;

figure(777)

while num<=Imax
    k=k+1;
    disp(100*num/Imax)
    
    for i=1:nz
        Cg1(:,:,i)=Cg1(:,:,i)/Wzg1(i);
    end
    Cg0(1:nx,1:ny,:)=Cg1;
    g0=real(sum(ifft2(Vg.*fft2(Cg0)),3));
    g00(1:nx,1:ny)=g1-g0(nx:(2*nx-1),ny:(2*ny-1));
    gf=fft2(g00); 
    
    for i=1:nz
        r1g0(:,:,i)=real(ifft2(Vg(:,:,i).*gf));
    end

    r1g=r1g0(nx:(2*nx-1),ny:(2*ny-1),:)-ag*Cg1;
    
    if k==1
        p1g(1:nx,1:ny,:)=r1g;
    else
        u2g=sum(r1g.*r1g,'all')/sum(r0g.*r0g,'all');
        p1g(1:nx,1:ny,:)=r1g+u2g*p1g(1:nx,1:ny,:);
    end
    r0g=r1g;
    %p0g=p1g; p0t=p1t;
    q1g=real(sum(ifft2(Vg.*fft2(p1g)),3)); 
    q11g=p1g(1:nx,1:ny,:);
    v2g=sum((r1g.*q11g),'all')/(sum(q1g(nx:(2*nx-1),ny:(2*ny-1)).*q1g(nx:(2*nx-1),ny:(2*ny-1)),'all')+ag*(sum(q11g.*q11g,'all')));
    Cg1=Cg1+v2g*q11g;
    for i=1:nz
        Cg1(:,:,i)=Cg1(:,:,i)*Wzg1(i);
    end
  
    
    rrg(k,1)=sum(r1g.*r1g,'all');
    plot(log10(rrg))
    num=num+1;
end
end