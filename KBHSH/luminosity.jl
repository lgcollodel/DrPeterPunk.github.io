include("radiant_energy.jl");


function get_oel2(xp,gm,rₕ,sqdetg,n)
    sp1 = Spline1D(xp,gm[:,1]); sp2 = Spline1D(xp,gm[:,2]); sp3 = Spline1D(xp,gm[:,3]);
    spg = Spline1D(xp,sqdetg)
    ϵ=0.0

    Ω    = zeros(n);
    Ωp   = zeros(n);
    Ωm   = zeros(n);
    Ωsp  = zeros(n);
    Ωsm  = zeros(n);
    en   = zeros(n);
    an   = zeros(n);
    a₊   = zeros(n);
    a₋   = zeros(n);
    v2p  = zeros(n);
    v2m  = zeros(n);
    Δr   = zeros(n);
    sqdg = zeros(n);
    sqdg[1]=sqdetg[1];
    tt=zeros(n);tp=zeros(n);pp=zeros(n);
    Ωp[1]=NaN; Ωm[1]=NaN;Ω[1]=NaN;v2p[1]=NaN;v2m[1]=NaN; Ωsp[1]=NaN; Ωsm[1]=NaN;
    xi=zeros(n);
    global Lflip = false;


    for i=2:n
        global Lvp  = false
        global Lvm  = false
        global Lup = false
        global Lum = false
        xi[i]=(i-1)/(n-1)
        x0 = xi[i]
        dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
        dx_dr2  = (((2*x0 + 1.)*(x0 - 1.)^2*rₕ^2 + 2*x0^3)*(x0 - 1.))/(sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*x0^2);
        gtt     = evaluate(sp1,x0);
        gtt_r   = derivative(sp1,x0)*dx_dr;
        gtt_rr  = derivative(sp1,x0,nu=2)*dx_dr^2+gtt_r*dx_dr2;
        gtp     = evaluate(sp2,x0);
        gtp_r   = derivative(sp2,x0)*dx_dr;
        gtp_rr  = derivative(sp2,x0,nu=2)*dx_dr^2+gtp_r*dx_dr2;
        gpp     = evaluate(sp3,x0);
        gpp_r   = derivative(sp3,x0)*dx_dr;
        gpp_rr  = derivative(sp3,x0,nu=2)*dx_dr^2+gpp_r*dx_dr2;
        Δ       = gtp_r^2-gtt_r*gpp_r;
        Δr[i]   = Δ
        sqdg[i] = evaluate(spg,x0);
        tt[i]=gtt;tp[i]=gtp;pp[i]=gpp;
        if Δ>=0
            Ωp[i] = (-gtp_r+sqrt(Δ))/gpp_r;
            Ωm[i] = (-gtp_r-sqrt(Δ))/gpp_r;
            Γ1      = -gtt-2*Ωp[i]*gtp-Ωp[i]^2*gpp;
            Γ2      = -gtt-2*Ωm[i]*gtp-Ωm[i]^2*gpp;
            if Γ1>0
                enp = -(gtt+gtp*Ωp[i])/sqrt(Γ1);
                anp = (gtp+gpp*Ωp[i])/sqrt(Γ1);
                v_rrp = enp^2*gpp_rr+2*enp*anp*gtp_rr+anp^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                v2p[i] = v_rrp
                if enp >= 0.0
                    v_rrp <= ϵ ? (global Lvp = true) : (global Lup = true);
                end
                #137 < i < 148 ? println(i,"\t",v_rrp) : nothing
            else
                enp=NaN; anp=NaN; v_rrp=NaN;
            end
            if Γ2>0
                enm = -(gtt+gtp*Ωm[i])/sqrt(Γ2);
                anm = (gtp+gpp*Ωm[i])/sqrt(Γ2);
                v_rrm = enm^2*gpp_rr+2*enm*anm*gtp_rr+anm^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                v2m[i] = v_rrm

                if enm >= 0.0
                    (Ωm[i] > 0.0 && v_rrm <= ϵ ) ? (global Lvm = true ) : (global Lum = true);
                end
                #118 < i < 138 ? println(i,"\t",v_rrm) : nothing
            else
                enm=NaN; anm=NaN; v_rrm=NaN;
            end
        else
            Ωp[i]=NaN; Ωm[i]=NaN; enp=NaN; enm=NaN; anp=NaN; anm=NaN; v_rrp=NaN; v_rrm=NaN;
        end
        a₊[i]=anp; a₋[i]=anm;
        if (Lvp && !Lvm)
            Ω[i]=Ωp[i]; en[i]=enp; an[i]=anp;
        elseif (Lvm && !Lvp)
            Ω[i]=Ωm[i]; en[i]=enm; an[i]=anm;
            println(i, "\t",xi[i],"\t", Ωp[i], "\t", Ωm[i])
        elseif (!Lvm && !Lvp)
            Ω[i]=NaN; en[i]=NaN; an[i]=NaN;
        elseif (Lvm && Lvp)
            Ω[i]=Ωp[i]; en[i]=enp; an[i]=anp;
        end
        Lvp ? Ωsp[i] = Ωp[i] : Ωsp[i] = NaN
        Lvm ? Ωsm[i] = Ωm[i] : Ωsm[i] = NaN
    end

    return xi, Ω,  en, an, sqdg, tt, tp, pp

end




function lum(isol,rh,mass,spin,ν,γ)
    #=
    gn = 6.6743015e-11
    cl = 2.99792458e8
    pc = 6.62607004e-34*gn/cl^3
    kᵦ = 1.38064852e-23*gn/cl^4
    σ  = 5.670374419e-8*gn/cl^5
    =#
    M₀  = 1.98847e33;
    yr  = 365*24*60*60.0;
    gn = 6.6743015e-8
    cl = 2.99792458e10
    pc = 6.62607004e-27
    kᵦ = 1.38064852e-16
    σ  = 5.670374419e-5
    α=2.5e6; β=2.e-6;
    factor=β*cl^6/(α*gn)^2/M₀/yr
    nr = 1001
    nϕ = 201
    if rh != 0
        xp,gd,sqd = get_metric(isol,rh);
        xi,o,en,an,sq,tt,tp,pp = get_oel2(xp,gd,rh,sqd,nr);
        is,ie = bound(o,xi);
        re = radiant_energy(is,ie,rh,xi,sq,en,an,o);
        fs = re*mass^2;
        rs  = r.(xi[is:ie],rh)/mass;
        lr = length(rs)
        o=o[is:ie];tt=tt[is:ie];tp=tp[is:ie];pp=pp[is:ie];
    elseif rh == 0
        rs, fs = nKerr_radiant_energy(spin,nr)
#        rs = riso.(rk,1,spin);
        lr = length(rs)
        tt = zeros(lr); tp = zeros(lr); pp=zeros(lr);
        o = 1.0./(spin .+ rs.^1.5)

        for i in 1:lr
            r0 = rs[i]
#            gk = kerr_metric(r0,M,spin)
            gk = nkerr_metric(r0,spin)
            tt[i] = gk[1]; tp[i] = gk[2]; pp[i] = gk[3]
        end
        M = α*M₀*gn/cl^2;
        #o = o*M;
        #rs=rs*M;
        fs=fs*factor;
    end
    temp = (fs/σ).^0.25
    println("T💕=",maximum(temp))
    ϕ = zeros(nϕ);
    redshift = zeros(lr,nϕ);
    intr = zeros(lr,nϕ);
    intϕ = zeros(nϕ);
    for j in 1:nϕ
        ϕ[j] = 2*π*(j-1)/(nϕ-1)
        for i in 1:lr
            redshift[i,j] = (1.0+o[i]*rs[i]*sin(γ)*sin(ϕ[j]))/sqrt(-tt[i]-2*o[i]*tp[i]-o[i]^2*pp[i])
#            redshift[i,j]=1.0
            intr[i,j] = ν^3*redshift[i,j]^3*rs[i]*M/(exp(ν*redshift[i,j]*pc/(kᵦ*temp[i]))-1.0)
        end
#        spr = Spline1D(rs,intr[:,j]/mass^2)
        spr = Spline1D(rs*M,intr[:,j]/mass^2)
        intϕ[j] = integrate(spr,rs[1]*M,rs[end]*M)
    end
    println("Integral in radial direction calculated")
    spϕ = Spline1D(ϕ, intϕ)
    luminosity = 8.0*cos(γ)/π*integrate(spϕ,ϕ[1],ϕ[end])/cl^2*pc
    println("ν=",ν,"\t","L=",luminosity,"\n")
    return luminosity
end
# fixed r (multiply by Mg)
# fixed factor
