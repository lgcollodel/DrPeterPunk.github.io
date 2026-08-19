include("sDifferenceFormula.jl")
include("sBS_funtions.jl");
# Declare Variables for the Run

Nd=1;  Nk=2;
Px=4;
M=[2,2];
Aleft=2 ; Aright=2;
Nlr=[Aleft,Aright];
# Get Initial guess
using DelimitedFiles, LinearAlgebra
#f0="h=0.102262840_ws=0.284000000000000_sigma=0.100";
f0="solitons_sol/h0=0.10159968466101595_om=0.2247_sigma=0.100";
#f0="phi=1.000000000000000_ws=0.812329116238928";
s0=readdlm(f0,Float64);
x=s0[:,1]; u0=s0[:,2:3];#u0=s0[:,4:5]
#p = [1.0,1.0,0.284,0.1]
p = [1.0,1.0,0.2247,0.1]
#p = [1.0,1.0,0.812329116238928,0.0]
potential=pot_solitonic;
#potential=pot_free;

#MAIN FUNCTION STARTS HERE
function main(Nd, Nk, Px, x, u0, p; r_tol=1.e-8)
    i0=0
    on=1.0
    u = SolutionVector(u0,u0);
    dform=DifferenceFormula(x,Px);
    arrsize=ArraySize(Nd,Nk,[dform.Nx]);
    du0=DerivArrays1D(arrsize); M=DerivMatrices1D(arrsize);
    du1=DerivArrays1D(arrsize);
    dform(M.dx,M.dxx);
    eq0=EqArray1D(arrsize);
    eq1=EqArray1D(arrsize);
    j=JacMatrices1D(arrsize);
    v0=reshape(u.u0,:,1);
    for k=1:Nk
        ux=@view du0.ux[:,k];uxx=@view du0.uxx[:,k]
        dform(ux,uxx,u.u0[:,k])
    end
    cjac!(j,u.u0,du0.ux,du0.uxx,M,p,x)
    for i0=1:dform.Nx
        equations(eq0.Pu,u.u0[i0,:],du0.ux[i0,:],du0.uxx[i0,:],p,x[i0])
        eq0.Res[i0] = eq0.Pu[1]
        eq0.Res[dform.Nx+i0]=eq0.Pu[2]
    end
    bc!(eq0.Res,u.u0,du0.ux)
    rmax=findmax(abs.(eq0.Res))[1]
    rmax < r_tol ? (println("Guess is the solution!"); return u0) : nothing
    while rmax > r_tol
        i0+=1
        if i0==1
            println(i0,"st Newtonian Iteration")
        elseif i0==2
            println(i0,"nd Newtonian Iteration")
        elseif i0==3
            println(i0,"rd Newtonian Iteration")
        else
            println(i0,"th Newtonian Iteration")
        end


        vn = - on*j.Jac\eq0.Res + v0

        u.u1=reshape(vn,(dform.Nx,2))
        for k=1:Nk
            sax=@view du1.ux[:,k];saxx=@view du1.uxx[:,k]
            dform(sax,saxx,u.u1[:,k])
        end
        for i0=1:dform.Nx
            equations(eq1.Pu,u.u1[i0,:],du1.ux[i0,:],du1.uxx[i0,:],p,x[i0])
            eq1.Res[i0] = eq1.Pu[1]
            eq1.Res[dform.Nx+i0]=eq1.Pu[2]
        end
        bc!(eq1.Res,u.u1,du1.ux)
        rmaxn=findmax(abs.(eq1.Res))[1]
        println("Err=",rmaxn)
        if rmaxn > rmax
            on=on/2
            on < 1.e-3 ? (println("No convergence");return 1) : nothing
        elseif (rmaxn < rmax) && (rmaxn > r_tol)

            rmax=rmaxn
            v0=vn
            eq0.Res=eq1.Res
            cjac!(j,u.u1,du1.ux,du1.uxx,M,p,x)
        elseif (rmaxn < r_tol)
            u.u0=u.u1;
            du0=du1;
            println("Convergence")
            return u.u0
        end
    end
end


for i=1:200
    println("Run n.",i,"\n ωₛ=",p[3])
    us = main(Nd,Nk,Px,x,u0,p)
    if us == 1
        break
    else
        u0 = us
        h0=u0[1,2]
        om=p[3]
        p[3]=p[3]-5.e-6
        println("h0=",h0)
        if i%100 == 0
            dform=DifferenceFormula(x,Px);
            ux=zeros(dform.Nx); uxx=zeros(dform.Nx);
            dform(ux,uxx,u0[:,1])
            M=ux[end]/2
            outf=string("h0=",h0,"_om=",om,"_sigma=0.100")
            io=open(outf,"a")
            for i=1:length(x)
                println(io,x[i],"\t",us[i,1],"\t",us[i,2])
            end
            close(io)
            println("M=",M)
        end
    end
end
function bc!(res,un,unx)
    Nx = length(un[:,1])
    res[1]=unx[1,1]
    res[Nx+1]=unx[1,2]

    res[Nx]=un[end,1]
    res[2*Nx]=un[end,2]
end

function cjac!(j,u0,ux,uxx,M,p,x)
    Nx   = length(x)
    Nk   = length(u0[1,:])
    Ju   = j.Ju
    Jux  = j.Jux
    Juxx = j.Juxx
    jac  = j.Jac
    Mx   = M.dx
    Mxx  = M.dxx
    for i0=1:Nx
        Du=@view Ju[i0,:,:]; Dux=@view Jux[i0,:,:]; Duxx=@view Juxx[i0,:,:];
        jacs(Du,Dux,Duxx,u0[i0,:],ux[i0,:],uxx[i0,:],p,x[i0])
    end
    for ke=1:Nk
        for kv=1:Nk
            jac[Nx*(ke-1)+1:ke*Nx,Nx*(kv-1)+1:kv*Nx] =
            Diagonal(Ju[1:Nx,ke,kv]) +
            Diagonal(Jux[1:Nx,ke,kv])*Mx +
            Diagonal(Juxx[1:Nx,ke,kv])*Mxx
        end
    end
    jac[1,:].=0.0; jac[1,1:Nx].=Mx[1,:];
    jac[(Nk-1)*Nx+1,:].=0.0; jac[(Nk-1)*Nx+1,(Nk-1)*Nx+1:Nk*Nx].=Mx[1,:];
    jac[Nx,:].=0.0; jac[Nx,Nx]=1.0;
    jac[2*Nx,:].=0.0; jac[2*Nx,2*Nx]=1.0;
end
