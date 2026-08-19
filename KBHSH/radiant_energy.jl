using DelimitedFiles, Dierckx, Glob, Plots;
include("get_metric.jl")



function rho(r0,M,J)
    (M+2*r0+J/M)*(M+2*r0-J/M)/(4*r0)
end
function riso(r0,M,a)
    (r0-M+sqrt(a^2-2*M*r0+r0^2))/2
end





function radiant_energy(is,ie,rₕ,x,sqdg,E,L,Ω)
    xs  = x[is];
    spo = Spline1D(x[is:ie],Ω[is:ie]);
    spl = Spline1D(x[is:ie],L[is:ie]);

    global Ω_r = zeros(ie-is+1);
    global L_r = zeros(ie-is+1);
    global integrand = zeros(ie-is+1);
    global RE = zeros(ie-is+1);

    for i=is:ie
        x0 = x[i]
        dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
        Ω_r[1+i-is] = derivative(spo,x0)*dx_dr
        L_r[1+i-is] = derivative(spl,x0)*dx_dr
        integrand[1+i-is] = (E[i]-L[i]*Ω[i])*L_r[1+i-is]/dx_dr
    end
    spi = Spline1D(x[is:ie],integrand);
    for i=is:ie
        (Ω_r[1+i-is] > 0 && abs(Ω_r[1+i-is]) < 1.e5) ? Ω_r[1+i-is] = 0.0 : nothing
         factor = -1.0/(4*π*sqdg[i])*Ω_r[1+i-is]/(E[i]-Ω[i]*L[i])^2
        RE[1+i-is] = factor*integrate(spi,xs,x[i])
    end
    return RE
end


function bound(Ω,xi)
    for i=1:length(Ω)
         if !isnan(Ω[i])
            global is=i
            break
        end
    end
    for i=length(xi):-1:1
        if xi[i]<0.97
            global ie=i
            break
         end
     end
     return is,ie
end





function Kerr_radiant_energy(M,a)
    # find isco
    n = 1001
    for i = 2:n
        rh = M+sqrt(M^2 - a^2)
        r0 = rh*(1.0 + (i-1)/(n-1)*24.0)
        gk = kerr_metric(r0,M,a)
        om = (a*M-r0^2*sqrt(M/r0))/(a^2*M-r0^3)
        d2 = -gk[1]-2*gk[2]*om-gk[3]*om^2
        if d2 > 0.0
            en = - (gk[1] + gk[2]*om)/sqrt(d2)
            an =   (gk[2] + gk[3]*om)/sqrt(d2)
            v_rr = en^2*gk[9]+2*an*en*gk[8]+an^2*gk[7] - (2*gk[8]*gk[2]+2*gk[5]^2-2*gk[4]*gk[6]-gk[7]*gk[3]-gk[1]*gk[9])
            v_rr <= 0.0 ? (global is = i; break) : nothing
        end
    end

    # get profile
    rc = zeros(n); Ωk = zeros(n); E = zeros(n); L = zeros(n);
    La = zeros(n); Ωa = zeros(n);
    L_r = zeros(n); Ωk_r = zeros(n); intk = zeros(n);
    for i = 1:n
        rh = M+sqrt(M^2 - a^2)
        r0 = rh*(1.0 + (i-2+is)/(n-1)*24.0)
        rc[i] = r0
        gk = kerr_metric(r0,M,a)
        Ωk[i] = (-a*M+sqrt(M*r0^3))/(r0^3-a^2*M)
        Ωa[i] = sqrt(M)/(a*sqrt(M)+r0^1.5)
#        Ωk_r[i] = 1.5*sqrt(M*r0)/(r0^3-a^2*M) - 3*r0^2*(-a*M+sqrt(M*r0^3))/(r0^3-a^2*M)^2
        Ωk_r[i] = -1.5*sqrt(M*r0)/(r0^1.5+a*sqrt(M))^2
        ds = sqrt(-gk[1]-2*gk[2]*Ωk[i]-gk[3]*Ωk[i]^2)
        E[i] = - (gk[1] + gk[2]*Ωk[i])/ds
        L[i] =   (gk[2] + gk[3]*Ωk[i])/ds
        La[i] =  sqrt(M)*(r0^2+a^2-2*a*sqrt(M*r0))/sqrt(r0^3-3*M*r0^2+2*a*r0*sqrt(M*r0))
        L_r[i] = (gk[5] + gk[6]*Ωk[i] + gk[3]*Ωk_r[i])/ds -
        L[i]/ds^2*(-gk[4] - 2*gk[5]*Ωk[i] - 2*gk[2]*Ωk_r[i] -
        gk[6]*Ωk[i]^2 - 2*gk[3]*Ωk[i]*Ωk_r[i])/2
        intk[i] = (E[i]-L[i]*Ωk[i])*L_r[i]
    end
    spo = Spline1D(rc,Ωk); spl = Spline1D(rc, L);
    # calculate integrand
#=
    for i = 1:n
        r0 = rc[i]
        Ωk_r[i] = derivative(spo,r0)
        L_r[i] = derivative(spl,r0)
        intk[i] = (E[i]-L[i]*Ωk[i])*L_r[i]
    end
=#
    spi = Spline1D(rc,intk)
    RK = zeros(n)
    for i =1:n
        factork = -Ωk_r[i]/(4*π*rc[i]*(E[i]-Ωk[i]*L[i])^2)
        RK[i] = factork*integrate(spi,rc[1],rc[i])
    end
    return rc/M, RK
end

function kerr_metric(r,M,a)
    gtt = -1+2*M/r; gtt_r = -2*M/r^2; gtt_rr = 4*M/r^3;
    gtp = -2*M*a/r; gtp_r = 2*a*M/r^2; gtp_rr = -4*a*M/r^3;
    gpp = r^2 + 2*M*a^2/r + a^2; gpp_r = 2*(r - a^2*M/r^2); gpp_rr = 2 +4*a^2*M/r^3;
    return [gtt,gtp,gpp,gtt_r,gtp_r,gpp_r,gtt_rr,gtp_rr,gpp_rr]
end


function nkerr_metric(r,j)
    gtt = -1+2/r; gtt_r = -2/r^2; gtt_rr = 4/r^3;
    gtp = -2*j/r; gtp_r = 2*j/r^2; gtp_rr = -4*j/r^3;
    gpp = r^2 + 2*j^2/r + j^2; gpp_r = 2*(r - j^2/r^2); gpp_rr = 2 +4*j^2/r^3;
    return [gtt,gtp,gpp,gtt_r,gtp_r,gpp_r,gtt_rr,gtp_rr,gpp_rr]
end
function nKerr_radiant_energy(j,n)
    # find isco
    #än = 1001
    for i = 2:n
        rh = 1.0+sqrt(1.0 - j^2)
        r0 = rh*(1.0 + (i-1)/(n-1)*24.0)
        gk = nkerr_metric(r0,j)
        om = 1.0/(j+r0^1.5)
        d2 = -gk[1]-2*gk[2]*om-gk[3]*om^2
        if d2 > 0.0
            en = - (gk[1] + gk[2]*om)/sqrt(d2)
            an =   (gk[2] + gk[3]*om)/sqrt(d2)
            v_rr = en^2*gk[9]+2*an*en*gk[8]+an^2*gk[7] - (2*gk[8]*gk[2]+2*gk[5]^2-2*gk[4]*gk[6]-gk[7]*gk[3]-gk[1]*gk[9])
            v_rr <= 0.0 ? (global is = i; break) : nothing
        end
    end

    # get profile
    rc = zeros(n); Ωk = zeros(n); E = zeros(n); L = zeros(n);
    L_r = zeros(n); Ωk_r = zeros(n); intk = zeros(n);
    for i = 1:n
        rh = 1.0+sqrt(1.0 - j^2)
        r0 = rh*(1.0 + (i-2+is)/(n-1)*24.0)
        rc[i] = r0
        gk = nkerr_metric(r0,j)
        Ωk[i] = 1.0/(j+r0^1.5)
        Ωk_r[i] = -1.5*sqrt(r0)/(j+r0^1.5)^2
        ds = sqrt(-gk[1]-2*gk[2]*Ωk[i]-gk[3]*Ωk[i]^2)
        E[i] = - (gk[1] + gk[2]*Ωk[i])/ds
        L[i] =   (gk[2] + gk[3]*Ωk[i])/ds
        L_r[i] = (gk[5] + gk[6]*Ωk[i] + gk[3]*Ωk_r[i])/ds -
        L[i]/ds^2*(-gk[4] - 2*gk[5]*Ωk[i] - 2*gk[2]*Ωk_r[i] -
        gk[6]*Ωk[i]^2 - 2*gk[3]*Ωk[i]*Ωk_r[i])/2
        intk[i] = (E[i]-L[i]*Ωk[i])*L_r[i]
    end
    spo = Spline1D(rc,Ωk); spl = Spline1D(rc, L);
    # calculate integrand
#=
    for i = 1:n
        r0 = rc[i]
        Ωk_r[i] = derivative(spo,r0)
        L_r[i] = derivative(spl,r0)
        intk[i] = (E[i]-L[i]*Ωk[i])*L_r[i]
    end
=#
    spi = Spline1D(rc,intk)
    RK = zeros(n)
    for i =1:n
        factork = -Ωk_r[i]/(4*π*rc[i]*(E[i]-Ωk[i]*L[i])^2)
        #factor = -1.0/(4*π*sqdg[i])*Ω_r[1+i-is]/(E[i]-Ω[i]*L[i])^2

        RK[i] = factork*integrate(spi,rc[1],rc[i])
    end
    return rc, RK
end



function Schw_Anal_rad(M,ṁ,r0)
    G = 6.6743015e-11
    c = 2.99792458e8
    mg = G/c^2*M

#    Ib = 1.0/sqrt(3*mg)*(log((sqrt(r0)-sqrt(3*mg))/(sqrt(r0)+sqrt(3*mg)))-log((sqrt(2)-1)/(sqrt(2)+1)))
#    Ia = sqrt(mg*r0)-sqrt(6)*mg-1.5*mg^1.5*Ib
    Ia = sqrt(mg)/2*(sqrt(3*mg)*log((sqrt(r0)-sqrt(3*mg))/(sqrt(r0)+sqrt(3*mg)))+2*sqrt(r0))
    F = 3*ṁ/(8*π*r0^1.5)*sqrt(mg)/(r0-3*mg)*Ia
    return F
end
