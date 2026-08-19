include("sDifferenceFormula.jl")

# Declare Variables for the Run

Nd=1;  Nk=2;
Px=4;
M=[2,2];
Aleft=2 ; Aright=2;
# Get Initial guess
using DelimitedFiles, LinearAlgebra
#f0="h=0.102262840_ws=0.284000000000000_sigma=0.100";
f0="h=0.10000000_ws=0.795753291469461_sigma=0.100"
s0=readdlm(f0,Float64);
x=s0[:,1]; u0=s0[:,4:5]
#p = [1.0,1.0,0.284,0.1]
p = [1.0,1.0,0.795753291469461,0.1]
function equations(Pu,u,ux,uxx,p,x)
    k1,xlam,om,sig = p
    f, h = u
    f_x,h_x = ux
    f_xx,h_xx = uxx
    r = x/(1.0-x)
    dxdr = (1.0-x)^2
    dxdr_x = 2*(x-1.0)


    Pu[1] =  - (exp(f)*(((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r
    ^2*sig^4*xlam - sig^4)*dxdr^3*f_x*h_x^2*r + 2*(h + sig)^2*(h - sig)^2*h^
    2*xlam)*k1*r - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*
    sig^4*xlam - sig^4)*dxdr_x*r - 2*sig^4)*dxdr*f_x - ((f_x^2 + f_xx)*(h^6*k1*
    r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4) + 2*
    (h + sig)^2*(h - sig)^2*h^2*h_x^2*k1^2*r^2*xlam)*dxdr^2*r) + ((f_xx + 4*
    h_x^2*k1)*dxdr^2*r^2 - 4 + (dxdr_x*r - 3)*dxdr*f_x*r)*h^2*k1*om^2*r*sig^4)


    Pu[2] = exp(f)*(((h^6*h_xx*k1*r^2*xlam - 3*h^5*h_x^2*k1*r^2*xlam - 2*h^4*
    h_xx*k1*r^2*sig^2*xlam + 4*h^3*h_x^2*k1*r^2*sig^2*xlam + h^2*h_xx*k1*r^2
    *sig^4*xlam - h*h_x^2*k1*r^2*sig^4*xlam - h_xx*sig^4 + (h^6*k1*r^2*xlam -
     2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*f_x*h_x)*dxdr^
    2 - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
    - sig^4)*dxdr^3*h_x^3*k1*r - (3*h^2 - sig^2)*(h + sig)*(h - sig)*h*xlam))*r
     + (2*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
     - sig^4)*h_x + (3*h^2 - sig^2)*(h + sig)*(h - sig)*f_x*h*r^2*xlam + (h^6*
    k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*
    dxdr_x*h_x*r)*dxdr) - ((f_x + h*h_x*k1 + dxdr_x*h*h_x*k1*r)*dxdr*r + (h*h_xx -
    h_x^2)*dxdr^2*k1*r^2 + 1)*h*om^2*r*sig^4
end
function jacs(Ju,Jux,Juxx,u,ux,uxx,p,x)
    k1,xlam,om,sig = p
    f, h = u
    f_x,h_x = ux
    f_xx,h_xx = uxx
    r = x/(1.0-x)
    dxdr = (1.0-x)^2
    dxdr_x = 2*(x-1.0)

    Ju[1,1] =  - exp(f)*(((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r
    ^2*sig^4*xlam - sig^4)*dxdr^3*f_x*h_x^2*r + 2*(h + sig)^2*(h - sig)^2*h^
    2*xlam)*k1*r - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*
    sig^4*xlam - sig^4)*dxdr_x*r - 2*sig^4)*dxdr*f_x - ((f_x^2 + f_xx)*(h^6*k1*
    r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4) + 2*
    (h + sig)^2*(h - sig)^2*h^2*h_x^2*k1^2*r^2*xlam)*dxdr^2*r)

    Jux[1,1] = (exp(f)*((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
    2*sig^4*xlam - sig^4)*dxdr_x*r - 2*sig^4 - (dxdr*h_x^2*k1*r - 2*f_x)*(h^6*
    k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*
    dxdr*r) - (dxdr_x*r - 3)*h^2*k1*om^2*r^2*sig^4)*dxdr

    Juxx[1,1] = (exp(f)*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
    2*sig^4*xlam - sig^4) - h^2*k1*om^2*r^2*sig^4)*dxdr^2*r

    Ju[1,2] =  - 2*(exp(f)*((dxdr^2*h_x^2*k1*r - dxdr_x)*dxdr*f_x*r^2 + 2 - (f_x^
    2 + f_xx + 2*h_x^2*k1)*dxdr^2*r^2)*(3*h^2 - sig^2)*(h + sig)*(h - sig)*xlam
     + ((f_xx + 4*h_x^2*k1)*dxdr^2*r^2 - 4 + (dxdr_x*r - 3)*dxdr*f_x*r)*om^2*sig
    ^4)*h*k1*r

    Jux[1,2] =  - 2*(exp(f)*((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*
    k1*r^2*sig^4*xlam - sig^4)*dxdr*f_x - 2*(h + sig)^2*(h - sig)^2*h^2*k1*r*
    xlam) + 4*h^2*k1*om^2*r*sig^4)*dxdr^2*h_x*k1*r^2

    Juxx[1,2] = 0

    Ju[2,1] = exp(f)*(((h^6*h_xx*k1*r^2*xlam - 3*h^5*h_x^2*k1*r^2*xlam - 2*h^4*
    h_xx*k1*r^2*sig^2*xlam + 4*h^3*h_x^2*k1*r^2*sig^2*xlam + h^2*h_xx*k1*r^2
    *sig^4*xlam - h*h_x^2*k1*r^2*sig^4*xlam - h_xx*sig^4 + (h^6*k1*r^2*xlam -
     2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*f_x*h_x)*dxdr^
    2 - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
    - sig^4)*dxdr^3*h_x^3*k1*r - (3*h^2 - sig^2)*(h + sig)*(h - sig)*h*xlam))*r
     + (2*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
     - sig^4)*h_x + (3*h^2 - sig^2)*(h + sig)*(h - sig)*f_x*h*r^2*xlam + (h^6*
    k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*
    dxdr_x*h_x*r)*dxdr)

    Jux[2,1] = (exp(f)*((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
    2*sig^4*xlam - sig^4)*dxdr*h_x + (3*h^2 - sig^2)*(h + sig)*(h - sig)*h*r*
    xlam) - h*om^2*r*sig^4)*dxdr*r

    Juxx[2,1] = 0

    Ju[2,2] = (exp(f)*((6*h^5*h_xx - 15*h^4*h_x^2 - 8*h^3*h_xx*sig^2 + 12*h^2*
    h_x^2*sig^2 + 2*h*h_xx*sig^4 - h_x^2*sig^4 + 2*(3*h^2 - sig^2)*(h + sig)*
    (h - sig)*f_x*h*h_x)*dxdr^2*k1*r^2 + 15*h^4 - 12*h^2*sig^2 + sig^4 - 2*(3*
    h^2 - sig^2)*(h + sig)*(h - sig)*dxdr^3*h*h_x^3*k1^2*r^3 + (2*(dxdr_x*r +
    2)*(3*h^2 - sig^2)*(h + sig)*(h - sig)*h*h_x*k1 + (15*h^4 - 12*h^2*sig^2 +
    sig^4)*f_x)*dxdr*r)*xlam - ((f_x + 2*h*h_x*k1 + 2*dxdr_x*h*h_x*k1*r)*dxdr*r + (
    2*h*h_xx - h_x^2)*dxdr^2*k1*r^2 + 1)*om^2*sig^4)*r

    Jux[2,2] = (exp(f)*((dxdr_x*r + 2 - 3*dxdr^2*h_x^2*k1*r^2)*(h^6*k1*r^2*xlam
    - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4) + ((h^6*k1*r
    ^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*f_x -
     2*(3*h^2 - sig^2)*(h + sig)*(h - sig)*h*h_x*k1*r^2*xlam)*dxdr*r) - ((dxdr_x*
    r + 1)*h - 2*dxdr*h_x*r)*h*k1*om^2*r^2*sig^4)*dxdr

    Juxx[2,2] = (exp(f)*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
    2*sig^4*xlam - sig^4) - h^2*k1*om^2*r^2*sig^4)*dxdr^2*r
end
function bc!(res,un,unx)
    res[1]=unx[1,1]
    res[dform.Nx+1]=unx[1,2]

    res[dform.Nx]=un[end,1]
    res[2*dform.Nx]=un[end,2]
end

#MAIN FUNCTION STARTS HERE
#function main(Nd, Nk, Px, x, u0, p)
# Incore:
dform=DifferenceFormula(x,Px);
arrsize=ArraySize(Nd,Nk,[dform.Nx]);



vu=DerivArrays1D(arrsize);
vun=DerivArrays1D(arrsize);




for k=1:Nk
    ux=@view vu.ux[:,k];uxx=@view vu.uxx[:,k]
    dform(ux,uxx,u0[:,k])
end




Mx=zeros(dform.Nx,dform.Nx); Mxx=zeros(dform.Nx,dform.Nx); dform(Mx,Mxx);
Pu=zeros(2); Ju=zeros(dform.Nx,2,2); Jux=zeros(dform.Nx,2,2); Juxx=zeros(dform.Nx,2,2);
dform(Mx,Mxx);




for i0=1:dform.Nx
    Du=@view Ju[i0,:,:]; Dux=@view Jux[i0,:,:]; Duxx=@view Juxx[i0,:,:];
    jacs(Du,Dux,Duxx,u0[i0,:],vu.ux[i0,:],vu.uxx[i0,:],p,x[i0])
end




jac=zeros(Nk*dform.Nx,Nk*dform.Nx);

for ke=1:Nk
    for kv=1:Nk
        jac[dform.Nx*(ke-1)+1:ke*dform.Nx,dform.Nx*(kv-1)+1:kv*(dform.Nx)] =
        Diagonal(Ju[1:dform.Nx,ke,kv]) +
        Diagonal(Jux[1:dform.Nx,ke,kv])*Mx +
        Diagonal(Juxx[1:dform.Nx,ke,kv])*Mxx
    end
end
#= Boundary Conditions
ux(0,1)=0; ux(0,3)=0;
u(end,1)=0; u(end,2)=0;
=#
jac[1,:].=0.0; jac[1,1:dform.Nx].=Mx[1,:];
jac[(Nk-1)*dform.Nx+1,:].=0.0; jac[(Nk-1)*dform.Nx+1,(Nk-1)*dform.Nx+1:Nk*dform.Nx].=Mx[1,:];
jac[dform.Nx,:].=0.0; jac[dform.Nx,dform.Nx]=1.0;
jac[2*dform.Nx,:].=0.0; jac[2*dform.Nx,2*dform.Nx]=1.0;


function cjac!(jac,Ju,Jux,Juxx,u0,ux,uxx,p,x)
    for i0=1:dform.Nx
        Du=@view Ju[i0,:,:]; Dux=@view Jux[i0,:,:]; Duxx=@view Juxx[i0,:,:];
        jacs(Du,Dux,Duxx,u0[i0,:],vu.ux[i0,:],vu.uxx[i0,:],p,x[i0])
    end
    for ke=1:Nk
        for kv=1:Nk
            jac[dform.Nx*(ke-1)+1:ke*dform.Nx,dform.Nx*(kv-1)+1:kv*(dform.Nx)] =
            Diagonal(Ju[1:dform.Nx,ke,kv]) +
            Diagonal(Jux[1:dform.Nx,ke,kv])*Mx +
            Diagonal(Juxx[1:dform.Nx,ke,kv])*Mxx
        end
    end
    jac[1,:].=0.0; jac[1,1:dform.Nx].=Mx[1,:];
    jac[(Nk-1)*dform.Nx+1,:].=0.0; jac[(Nk-1)*dform.Nx+1,(Nk-1)*dform.Nx+1:Nk*dform.Nx].=Mx[1,:];
    jac[dform.Nx,:].=0.0; jac[dform.Nx,dform.Nx]=1.0;
    jac[2*dform.Nx,:].=0.0; jac[2*dform.Nx,2*dform.Nx]=1.0;
end



for k=1:Nk
    sax=@view vu.ux[:,k];saxx=@view vu.uxx[:,k]
    dform(sax,saxx,u0[:,k])
end
cjac!(jac,Ju,Jux,Juxx,u0,vu.ux,vu.uxx,p,x)
res=zeros(2*dform.Nx);resn=zeros(2*dform.Nx);

for i0=1:dform.Nx
    equations(Pu,u0[i0,:],vu.ux[i0,:],vu.uxx[i0,:],p,x[i0])
    res[i0] = Pu[1]
    res[dform.Nx+i0]=Pu[2]
end
bc!(res,u0,vu.ux)
rmax=findmax(abs.(res))[1]
rtol=1.e-8
i0=0
on=1.0
un=zeros(dform.Nx);
v0=reshape(u0,:,1);
while rmax > rtol
    i0+=1
    vn = - on*jac\res + v0

    un=reshape(vn,(dform.Nx,2))
    for k=1:Nk
        sax=@view vun.ux[:,k];saxx=@view vun.uxx[:,k]
        dform(sax,saxx,un[:,k])
    end
    for i0=1:dform.Nx
        equations(Pu,un[i0,:],vun.ux[i0,:],vun.uxx[i0,:],p,x[i0])
        resn[i0] = Pu[1]
        resn[dform.Nx+i0]=Pu[2]
    end
    bc!(resn,un,vun.ux)
    rmaxn=findmax(abs.(resn))[1]

    println(i0)
    if rmaxn > rmax
        on=on/2
        on < 1.e-3 ? (println("No convergence");break) : nothing
    else
        rmax=rmaxn
        v0=vn
        res=resn
        cjac!(jac,Ju,Jux,Juxx,un,vun.ux,vun.uxx,p,x)
    end
end
if rmaxn < rtol
    u0=un;
    vu=vun;
    println("Convergence")
end
