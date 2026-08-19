include("sDifferenceFormula.jl")
include("functions.jl")
using LinearAlgebra
# Ex 1. u'' + 2 u' + 2 = 0
Nd=1;  Nk=1;
Px=2;
M=[2];

Nx=31;
x=collect(range(0,1,length=Nx));
dform=DifferenceFormula(x,Px);

Mx=zeros(dform.Nx,dform.Nx); Mxx=zeros(dform.Nx,dform.Nx);
dform(Mx,Mxx)

function eq(u,ux,uxx)
    return uxx+2*ux+u
end


u0 = 2 .-x;
u0x = Mx*u0; u0xx = Mxx*u0;

res = eq.(u0,u0x,u0xx)
res[1]=u0[1]-1
res[end]=u0[end]
using LinearAlgebra

jac = I + 2*Mx + Mxx;
# bc u(0)=1; u(1)=0;
jac[1,1]=1; jac[1,2:end].=0;
jac[end,end]=1; jac[end,1:end-1].=0;
b=zeros(Nx); b[1]=1;

rmax=findmax(abs.(res))[1]
rtol=1.e-5
i0=0
while rmax > rtol
    i0+=1
    un = - jac\res + u0

    resn= jac*un
    resn[1]=un[1]-1
    resn[end]=un[end]
    rmaxn=findmax(abs.(resn))[1]
    println(i0)
    if rmaxn > rmax
        println("DIVERGENCE")
        break
    else
        rmax=rmaxn
        u0=un
        res=resn
    end
end


###############################
# Ex 2. u'' + 2 x u'^3 = 0
Nd=1;  Nk=1;
Px=2;
M=[2];
Nx=31;
Nlr=[1,1];
x=collect(range(1,exp(1),length=Nx))
a=1/(exp(1)-1); b=-a;
u0=zeros(Nx,Nk);
u0[:,1] = a*x.+b;
p=[];

function ueqs!(Pu,u,ux,uxx,p,x)
    Pu[1] = uxx[1]+x*ux[1]^3
end
function ubc!(c,u,ux,p)
    c[1] = u[1]
    c[2] = u[1] - 1.0
end
function ujac_eq!(Ju,Jux,Juxx,u,ux,uxx,p,x)
    Ju[1,1] = 0.0
    Jux[1,1] = 3*ux[1]^2*x
    Juxx[1,1]=1.0
end
function ujac_bc!(dc,dcx,u,ux)
    dc[1,1] = 1.0
    dc[2,1] = 1.0
end

us = main(Nd,Px,x,u0,p,Nlr)



u0x = Mx*u0; u0xx = Mxx*u0;

res = eq.(u0,u0x,u0xx,x)
res[1]=u0[1];
res[end]=u0[end]-1.0;


jac = Mxx + Diagonal(3*u0x.^2)*Diagonal(x)*Mx
# bc u(x0)=0; u(x1)=1;
jac[1,1]=1; jac[1,2:end].=0;
jac[end,end]=1; jac[end,1:end-1].=0;


b=zeros(Nx); b[1]=1;

rmax=findmax(abs.(res))[1]
rtol=1.e-12
i0=0
while rmax > rtol
    i0+=1;
    un = - jac\res + u0;

    unx = Mx*un; unxx = Mxx*un;
    res1 = eq.(un,unx,unxx,x);
    res1[1]=u0[1];
    res1[end]=u0[end]-1.0;
    rmax1=findmax(abs.(res1))[1]
    resn= jac*un;
    resn[1]=un[1];
    resn[end]=un[end]-1;
    rmaxn=findmax(abs.(resn))[1]
    println(i0)
    if rmax1 > rmax
        println("DIVERGENCE")
        break
    else
        rmax=rmax1;
        u0=un;
        u0x = unx; u0xx = unxx;
        jac = Mxx + Diagonal(3*u0x.^2)*Diagonal(x)*Mx;
        jac[1,1]=1; jac[1,2:end].=0;
        jac[end,end]=1; jac[end,1:end-1].=0;
        res=res1;
    end
end
