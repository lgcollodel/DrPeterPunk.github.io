# FIX INDICES TO CORRECT RES[I] WHEN NONSYMETRIC CONDITIONS
# For multiple conditions for a single field the marix is still sparse but no longe block-banded
function boundary_conditions!(res,u,ux,Nlr,ulr,p)
    c=zeros(Nlr[1]+Nlr[2]);
    Nx=length(u[:,1])
    for i=1:Nlr[1]
    	j=(ulr[1][i]-1)*Nx+1
        ubc!(c,u[1,:],ux[1,:])
        res[j]=c[i]
    end
    for i=1:Nlr[2]#i=Nlr[1]+1:Nlr[2]
        ii=i+Nlr[1]
        j=ulr[2][i]*Nx
        ubc!(c,u[end,:],ux[end,:])
        res[j]=c[ii]
    end
end
function jac!(j,u,ux,uxx,M,p,x,Nlr,ulr)
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
        ujac_eq!(Du,Dux,Duxx,u0[i0,:],ux[i0,:],uxx[i0,:],p,x[i0])
    end
    for ke=1:Nk
        for kv=1:Nk
            jac[Nx*(ke-1)+1:ke*Nx,Nx*(kv-1)+1:kv*Nx] =
            Diagonal(Ju[1:Nx,ke,kv]) +
            Diagonal(Jux[1:Nx,ke,kv])*Mx +
            Diagonal(Juxx[1:Nx,ke,kv])*Mxx
        end
    end
    # Boundaries
    dc=zeros(Nlr[1]+Nlr[2],Nk);
    dcx=zeros(Nlr[1]+Nlr[2],Nk);
    ujac_bc!(dc,dcx,u,ux)
    for i=1:Nlr[1]
    	j=(ulr[1][i]-1)*Nx+1
        jac[j,:].=0.0
        for k=1:Nk
            ik=1+Nx*(k-1)
            jac[j,ik:Nx+ik-1].=dcx[i,k]*Mx[1,:];
            jac[j,ik]+=dc[i,k]
        end
    end
    for i=1:Nlr[2]
        ii=i+Nlr[1]
        j=ulr[2][i]*Nx
        jac[j,:].=0.0
        for k=1:Nk
            ik=1+Nx*(k-1)
            jac[j,ik:Nx+ik-1].=dcx[ii,k]*Mx[end,:];
            jac[j,Nx*k]+=dc[ii,k]
        end
    end
end
function main(Nd, Px, x, u0, p, Nlr,ulr; r_tol=1.e-8)
    i0=0
    on=1.0
    Nk = length(u0[1,:])
    Nx = length(u0[:,1])
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
    jac!(j,u.u0,du0.ux,du0.uxx,M,p,x,Nlr,ulr)
    for i0=1:dform.Nx
        ueqs!(eq0.Pu,u.u0[i0,:],du0.ux[i0,:],du0.uxx[i0,:],p,x[i0])
        for k=1:Nk
            ik=i0+Nx*(k-1)
            eq0.Res[ik] = eq0.Pu[k]
        end
        #eq0.Res[i0] = eq0.Pu[1]
        #eq0.Res[dform.Nx+i0]=eq0.Pu[2]
    end
    boundary_conditions!(eq0.Res,u0,du0.ux,Nlr,ulr,p)
    rmax=findmax(abs.(eq0.Res))[1]
    rmax < r_tol ? (println("Guess is the solution!"); return u0) : nothing
    while rmax > r_tol
        i0+=1

        vn = - on*j.Jac\eq0.Res + v0

        Nk > 1 ? u.u1=reshape(vn,(dform.Nx,Nk)) : u.u1 = vn
        for k=1:Nk
            sax=@view du1.ux[:,k];saxx=@view du1.uxx[:,k]
            dform(sax,saxx,u.u1[:,k])
        end
        for i0=1:dform.Nx
            ueqs!(eq1.Pu,u.u1[i0,:],du1.ux[i0,:],du1.uxx[i0,:],p,x[i0])
            for k=1:Nk
                ik=i0+Nx*(k-1)
                eq1.Res[ik] = eq1.Pu[k]
            end
            #eq1.Res[i0] = eq1.Pu[1]
            #eq1.Res[dform.Nx+i0]=eq1.Pu[2]
        end
        boundary_conditions!(eq1.Res,u.u1,du1.ux,Nlr,ulr,p)
        rmaxn=findmax(abs.(eq1.Res))[1]
        !isfinite(rmaxn) ? (println("Newton Method Diverges"); return 1) : nothing
        if rmaxn > rmax
            on=on/2
            on < 1.e-3 ? (println("No convergence in ",i0," iterations \n Maximum of the residual: ",rmaxn);return 1) : nothing
        elseif (rmaxn < rmax) && (rmaxn > r_tol)

            rmax=rmaxn
            v0=vn
            eq0.Res=eq1.Res
            jac!(j,u.u1,du1.ux,du1.uxx,M,p,x,Nlr,ulr)
        elseif (rmaxn < r_tol)
            u.u0=u.u1;
            du0=du1;
            #println("Convergence")
            return u.u0
        end
    end
end
