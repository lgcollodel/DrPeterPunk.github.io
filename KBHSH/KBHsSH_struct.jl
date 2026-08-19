#=
This general routine is fed with
In: [String] solution
In: [Float64] Horizon radius
In: [Float64] Natural frequency
In: [Int64] Winding Number (default m=1)
Note that because of construction, usually the BH rotates in the negative direction.
For this reason we flip it when reading, so the sign changes. This is also reflected in the
equations for the Komar expressions where change the sign of m. Be careful with the notation.
=#

using DelimitedFiles, Dierckx;

mutable struct KBHsSH <: Function
    solution :: String
    rₕ       :: Float64
    m        :: Int64
    ωₛ       :: Float64
    ADM_M    :: Float64
    ADM_J    :: Float64
    Hole_M   :: Float64
    Hole_J   :: Float64
    Hair_M   :: Float64
    Hair_J   :: Float64
    charge   :: Float64
    x        :: Array{Float64,1}
    r        :: Array{Float64,1}
    θ        :: Array{Float64,1}
    gm       :: Array{Float64,3}
    gm_eq    :: Array{Float64,2}
    sqdetg   :: Array{Float64,2}


    function KBHsSH(solution,rₕ,ωₛ;m=1)
        self = new()
        self.solution = solution
        self.rₕ = rₕ
        self.ωₛ = ωₛ
        self.m  = m
        metric_data = metric(solution,rₕ,ωₛ;m=m)
        self.r, self.x, self.θ, self.gm, self.sqdetg = metric_data[1:5]
        self.ADM_M, self.ADM_J, self.Hair_M, self.Hair_J, self.Hole_M, self.Hole_J = metric_data[6:11]
        self.gm_eq = @view self.gm[:,end,:]
        self.charge = self.Hair_J/(self.Hair_J+self.Hole_J)
        return self
    end
end

mutable struct Circular_Orbit <: Function

    prograde :: Bool

    r :: Array{Float64,1}
    Ω :: Array{Float64,1}
    E :: Array{Float64,1}
    L :: Array{Float64,1}

    circular_orbit :: Array{Bool,1}
    timelike       :: Array{Bool,1}
    stable         :: Array{Bool,1}

    function Circular_Orbit(r,gm_spline;prograde=true)
        self = new()
        Nr = length(r)
        self.r = r
        self.prograde = prograde

        self.Ω = Array{Float64,1}(undef,Nr)
        self.E = Array{Float64,1}(undef,Nr)
        self.L = Array{Float64,1}(undef,Nr)

        self.circular_orbit = Array{Bool,1}(undef,Nr)
        self.timelike       = Array{Bool,1}(undef,Nr)
        self.stable         = Array{Bool,1}(undef,Nr)

        params = [0.0,0.0,0.0,false,false,false]
        for i=1:Nr
            circular_orbit!(params, gm_spline,r[i]; prograde)
            self.Ω[i] = params[1]
            self.E[i] = params[2]
            self.L[i] = params[3]

            self.circular_orbit[i] = params[4]
            self.timelike[i]       = params[5]
            self.stable[i]         = params[6]
        end
        return self
    end
end


function x_to_r(x,rₕ)
    sqrt((x - 1.0)^2*rₕ^2 + x^2)/(1.0-x)
end

function r_to_x(r,rₕ)
    - ((r + rₕ)*(r - rₕ) - sqrt(r^2 - rₕ^2))/(rₕ^2 + 1.0 - r^2)
end
#=
function rk_to_r(rh,M,a,rk)
    dma=M^2-a^2

    return (rh*(sqrt(dma)*sqrt( - dma - 2*M*rk + M^2 + rk^2) - sqrt( - dma - 2*M*rk + M^2 +
    rk^2)*M + sqrt( - dma - 2*M*rk + M^2 + rk^2)*rk - sqrt(dma)*M + sqrt(dma)*rk -
    2*M*rk + M^2 + rk^2))/(2*sqrt(dma)*(sqrt( - dma - 2*M*rk + M^2 + rk^2) - M + rk))
end
=#
function rk_to_r(rh,M,a,rk)
    dma=M^2-a^2
    return (rh*( - sqrt(dma)*M + sqrt(dma)*rk + dma))/(2*dma)
end
function drk_dr(rh,M,a)
    return drkdr=2*sqrt(M^2-a^2)/rh
end
function r_to_rk(rh,M,a,r)
    dma=M^2-a^2
    return (sqrt(dma)*M*rh + 2*dma*r - dma*rh)/(sqrt(dma)*rh)
end
function kerr_to_r(M,a,r)
    return r-a^2/(M+sqrt(M^2-a^2))
end
function r_to_kerr(M,a,r)
    return r + a^2/(M+sqrt(M^2-a^2))
end
function mass_rh(j,rₕ)
    return rₕ/2*(1+sqrt(1-j^2))/(1-j^2+sqrt(1-j^2))
end

function metric(isol,rₕ,ωₛ;m=1)
    sol = readdlm(isol,Float64);
    xx = @view sol[:,1];
    Nxt = length(xx);
    x=Float64[]
    for i=1:Nxt
        push!(x,xx[i])
        xx[i]==1.0 ? break : nothing
    end
    r = x_to_r.(x,rₕ)
    Nx = length(x)
    Ny = Int(Nxt/Nx)
    i_eq = (Ny-1)*Nx+1;

    θ  = Array{Float64,1}(undef,Ny)
    gm     = Array{Float64,3}(undef,Nx,Ny,6)
    mfunc  = Array{Float64,3}(undef,Nx,Ny,4)
    sqdetg = Array{Float64,2}(undef,Nx,Ny)
    dmv    = Array{Float64,1}(undef,Nx)
    djv    = Array{Float64,1}(undef,Nx)
    dmh    = Array{Float64,1}(undef,Ny)
    djh    = Array{Float64,1}(undef,Ny)
    dmy    = Array{Float64,1}(undef,Ny)
    djy    = Array{Float64,1}(undef,Ny)

    for j=1:Ny
        θ[j] = sol[(j-1)*Nx+1,2]
        st = sin(θ[j])
        for i=1:Nx
            r0 = r[i]
            line =(j-1)*Nx+i
            mfunc[i,j,1] = f0 =  sol[line,3]
            mfunc[i,j,2] = f1 =  sol[line,4]
            mfunc[i,j,3] = f2 =  sol[line,5]
            mfunc[i,j,4] = ω  = -sol[line,6]


            e2f0 = exp(2*f0)
            e2f1 = exp(2*f1)
            e2f2 = exp(2*f2)

            gm[i,j,1] = -(1.0-rₕ/r0)*e2f0+e2f2*ω^2*st^2;
            gm[i,j,2] = e2f1/(1.0-rₕ/r0)
            gm[i,j,3] = e2f1*r0^2
            gm[i,j,4] = -e2f2*ω*r0*st^2
            gm[i,j,5] = e2f2*r0^2*st^2
            gm[i,j,6] = ϕ  =  sol[line,7]

            sqdetg[i,j] = exp(f0+f2)*e2f1*r0^2*st

            if (i==1 || i==Nx || j==1)
                dmv[i]=0.0
                djv[i]=0.0
            else
                dxdr   = r0*(1.0-x[i])^3/x[i]
                dmv[i] = -sqdetg[i,j]/e2f0*(-2*(-ω*m+ωₛ*r0)*ωₛ + (r0-rₕ)*e2f0)*ϕ^2/((r0-rₕ)*dxdr)
                djv[i] = sqdetg[i,j]/e2f0*(-ω*m+ωₛ*r0)*ϕ^2/((r0-rₕ)*dxdr)
            end
        end
        spo = Spline1D(x,mfunc[:,j,4])
        ωxx = derivative(spo,0,nu=2)
        f0  = mfunc[1,j,1]
        f2  = mfunc[1,j,3]
        ω   = mfunc[1,j,4]
        dmh[j] = 0.5*exp(f2-f0)*((ω-ωxx*rₕ^2)*exp(2*f2)*ω*st^2+exp(2*f0))*rₕ*st
        djh[j] = 0.25*exp(f2-f0)*(ω-ωxx*rₕ^2)*exp(2*f2)*rₕ^2*st^3


        spm = Spline1D(x,dmv)
        dmy[j] = integrate(spm,0,1)


        spj = Spline1D(x,djv)
        djy[j] = integrate(spj,0,1)

    end


    smv = Spline1D(θ, dmy)
    sjv = Spline1D(θ, djy)
    smh = Spline1D(θ, dmh)
    sjh = Spline1D(θ, djh)

    mv =  integrate(smv, 0, π/2)
    jv =  integrate(sjv, 0, π/2)
    mh =  integrate(smh, 0, π/2)
    jh =  integrate(sjh, 0, π/2)

    f0_eq = @view mfunc[:,end,1]
    ω_eq  = @view mfunc[:,end,4]


    sf0 = Spline1D(x,f0_eq)
    sω  = Spline1D(x,ω_eq)
    mass = derivative(sf0,1.0) + rₕ/2
    angm = derivative(sω,1.0,nu=2)/4

    metric_data = [r, x, θ, gm, sqdetg, mass, angm, mv, jv, mh, jh]
    return metric_data
end

function metric2D(isol,rₕ,ωₛ;m=1)

        sol = readdlm(isol,Float64);
        xx = @view sol[:,1];
        Nxt = length(xx);
        x=Float64[]
        for i=1:Nxt
            push!(x,xx[i])
            xx[i]==1.0 ? break : nothing
        end
        r = x_to_r.(x,rₕ)
        Nx = length(x)
        Ny = Int(Nxt/Nx)
        i_eq = (Ny-1)*Nx+1;

        θ  = Array{Float64,1}(undef,Ny)
        gm     = Array{Float64,3}(undef,Nx,Ny,6)
        mfunc  = Array{Float64,3}(undef,Nx,Ny,4)
        sqdetg = Array{Float64,2}(undef,Nx,Ny)
        dmv    = Array{Float64,2}(undef,Nx,Ny)
        djv    = Array{Float64,2}(undef,Nx,Ny)
        dmh    = Array{Float64,1}(undef,Ny)
        djh    = Array{Float64,1}(undef,Ny)

        for j=1:Ny
            θ[j] = sol[(j-1)*Nx+1,2]
            st = sin(θ[j])
            for i=1:Nx
                r0 = r[i]
                line =(j-1)*Nx+i
                mfunc[i,j,1] = f0 =  sol[line,3]
                mfunc[i,j,2] = f1 =  sol[line,4]
                mfunc[i,j,3] = f2 =  sol[line,5]
                mfunc[i,j,4] = ω  = -sol[line,6]


                e2f0 = exp(2*f0)
                e2f1 = exp(2*f1)
                e2f2 = exp(2*f2)

                gm[i,j,1] = -(1.0-rₕ/r0)*e2f0+e2f2*ω^2*st^2;
                gm[i,j,2] = e2f1/(1.0-rₕ/r0)
                gm[i,j,3] = e2f1*r0^2
                gm[i,j,4] = -e2f2*ω*r0*st^2
                gm[i,j,5] = e2f2*r0^2*st^2
                gm[i,j,6] = ϕ  =  sol[line,7]

                sqdetg[i,j] = exp(f0+f2)*e2f1*r0^2*st

                if (i==1 || i==Nx || j==1)
                    dmv[i]=0.0
                    djv[i]=0.0
                else
                    dxdr   = r0*(1.0-x[i])^3/x[i]
                    dmv[i,j] = -sqdetg[i,j]/e2f0*(-2*(-ω*m+ωₛ*r0)*ωₛ + (r0-rₕ)*e2f0)*ϕ^2/((r0-rₕ)*dxdr)
                    djv[i,j] = sqdetg[i,j]/e2f0*(-ω*m+ωₛ*r0)*ϕ^2/((r0-rₕ)*dxdr)
                end
            end
            spo = Spline1D(x,mfunc[:,j,4])
            ωxx = derivative(spo,0,nu=2)
            f0  = mfunc[1,j,1]
            f2  = mfunc[1,j,3]
            ω   = mfunc[1,j,4]
            dmh[j] = 0.5*exp(f2-f0)*((ω-ωxx*rₕ^2)*exp(2*f2)*ω*st^2+exp(2*f0))*rₕ*st
            djh[j] = 0.25*exp(f2-f0)*(ω-ωxx*rₕ^2)*exp(2*f2)*rₕ^2*st^3

        end


        smv = Spline2D(x,θ,dmv)
        sjv = Spline2D(x,θ,djv)


        smh = Spline1D(θ, dmh)
        sjh = Spline1D(θ, djh)

        mv =  integrate(smv, 0, 1, 0, π/2)
        jv =  integrate(sjv, 0, 1, 0, π/2)
        mh =  integrate(smh, 0, π/2)
        jh =  integrate(sjh, 0, π/2)

        f0_eq = @view mfunc[:,end,1]
        ω_eq  = @view mfunc[:,end,4]


        sf0 = Spline1D(x,f0_eq)
        sω  = Spline1D(x,ω_eq)
        mass = derivative(sf0,1.0) + rₕ/2
        angm = derivative(sω,1.0,nu=2)/4

        metric_data = [r, x, θ, gm, sqdetg, mass, angm, mv, jv, mh, jh]
        return metric_data
end

function circular_orbit!(params, gm_spline,r; prograde=true)
    Ω, E, L, circular_orbit, timelike, stable = params
    gtt = evaluate(gm_spline[1],r); gtt_r = derivative(gm_spline[1],r); gtt_rr = derivative(gm_spline[1],r,nu=2);
    gtp = evaluate(gm_spline[2],r); gtp_r = derivative(gm_spline[2],r); gtp_rr = derivative(gm_spline[2],r,nu=2);
    gpp = evaluate(gm_spline[3],r); gpp_r = derivative(gm_spline[3],r); gpp_rr = derivative(gm_spline[3],r,nu=2);

    Δ =  - gpp_r*gtt_r + gtp_r^2
    Δ < 0 ? circular_orbit = false : circular_orbit = true

    if circular_orbit
        if prograde
            Ω = (sqrt(Δ) - gtp_r)/gpp_r
        else
            Ω = -(sqrt(Δ) + gtp_r)/gpp_r
        end
    else
        Ω = NaN
    end

    Γ=gpp*Ω^2 + 2*gtp*Ω + gtt
    Γ > 0 ? timelike = false : timelike = true

    if timelike
        E = - sqrt(( - 1)/Γ)*(gtp*Ω + gtt)
        L =   sqrt(( - 1)/Γ)*(gpp*Ω + gtp)
    else
        E = NaN; L = NaN;
    end

    #  Check Stability
    if circular_orbit && timelike
        v_rr = (E^2*gpp_rr+2*E*L*gtp_rr+L^2*gtt_rr
        -2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr);
        v_rr >=0.0 ? stable = false : stable = true
    else
        stable = false
    end

    params[:] = [Ω, E, L, circular_orbit, timelike, stable]
end

function effective_potential!(V,gm,L; ϵ=1)
    Nr = length(V)
    for i = 2:Nr-1
        Δ = gm[i,4]^2 - gm[i,1]*gm[i,5]
        A = -L*gm[i,4]
        B = sqrt(Δ*(gm[i,5]+L^2))
        V[i] = (A + ϵ*B)/gm[i,5]
    end
    V[1]= -L*gm[1,4]/gm[1,5]
    V[end] = 1.0
end


function r̈(params,gm_spline,r)
    E, L = params
    gtt = evaluate(gm_spline[1],r); gtt_r = derivative(gm_spline[1],r);
    gtp = evaluate(gm_spline[2],r); gtp_r = derivative(gm_spline[2],r);
    gpp = evaluate(gm_spline[3],r); gpp_r = derivative(gm_spline[3],r);


    return ((((gtp*gtt_r - 2*gtp_r*gtt)*gtp + gpp_r*gtt^2)*grr - Δ*grr_r*gtt)*L^2 + (
    Δ - E^2*gpp)*Δ*grr_r + (gpp^2*gtt_r - 2*gpp*gtp*gtp_r + gpp_r*gtp^2)*E
    ^2*grr + 2*((gpp_r*gtt - 2*gtp*gtp_r + gpp*gtt_r)*grr*gtp + (grr*gtp_r - grr_r*
    gtp)*Δ)*L*E)/(2*Δ^2*grr^2)

end

function orbit!(du,u,p,t)
    E = p[1]; L = p[2];
    gm_spline = p[3:6];
    r = u[3]
    gtt = evaluate(gm_spline[1],r); gtt_r = derivative(gm_spline[1],r);
    grr = evaluate(gm_spline[2],r); grr_r = derivative(gm_spline[2],r);
    gtp = evaluate(gm_spline[3],r); gtp_r = derivative(gm_spline[3],r);
    gpp = evaluate(gm_spline[4],r); gpp_r = derivative(gm_spline[4],r);
    Δ = gtp^2 - gtt*gpp
    du[1] = (L*gtp + E*gpp)/Δ
    du[2] = -(L*gtt + E*gtp)/Δ
    du[3] = u[4]

    du[4] = ((((gtp*gtt_r - 2*gtp_r*gtt)*gtp + gpp_r*gtt^2)*grr - Δ*grr_r*gtt)*L^2 + (
    Δ - E^2*gpp)*Δ*grr_r + (gpp^2*gtt_r - 2*gpp*gtp*gtp_r + gpp_r*gtp^2)*E
    ^2*grr + 2*((gpp_r*gtt - 2*gtp*gtp_r + gpp*gtt_r)*grr*gtp + (grr*gtp_r - grr_r*
    gtp)*Δ)*L*E)/(2*Δ^2*grr^2)

end


function ic!(u,params,gm_spline,r)
    E, L = params
    gtt = evaluate(gm_spline[1],r);
    grr = evaluate(gm_spline[2],r);
    gtp = evaluate(gm_spline[3],r);
    gpp = evaluate(gm_spline[4],r);
    Δ = gtp^2 - gtt*gpp

    u[1] = 0.0
    u[2] = 0.0
    u[3] = r
    u[4] = sqrt(((L*gtt + 2*E*gtp)*L - (Δ - E^2*gpp))/(Δ*grr))
end




#=
function line!(ab::Array{Float64,1},x1,x2,y1,y2)
       a=(y2-y1)/(x2-x1)
       b=y1-a*x1
       ab[1]=a
       ab[2]=b
end
ab=[0.0,0.0]
rst=0.933442507211974
#i=59
istart=2
iend=2
lr=length(black_hole.r)-iend
rb=black_hole.r[istart:lr]
rk=r_to_rk.(rₕ,black_hole.ADM_M,black_hole.ADM_J/black_hole.ADM_M,black_hole.r);
rhk=r_to_rk.(rₕ,black_hole.ADM_M,black_hole.ADM_J/black_hole.ADM_M,rₕ)
#rk=r_to_rk.(rₕ,black_hole.ADM_M,black_hole.ADM_J/black_hole.ADM_M,black_hole.r[istart:lr]);
#line!(ab,rb[1],rb[50],rk[1],rk[50])
#α=ab[1]
#rk2=ab[1]*rb.+ab[2]
#r0=rk[55]
#=
sptt = Spline1D(rk,black_hole.gm_eq[istart:lr,1]); # Needs to discard infinity
sprr = Spline1D(rk,black_hole.gm_eq[istart:lr,2]);
sptφ = Spline1D(rk,black_hole.gm_eq[istart:lr,4]);
spφφ = Spline1D(rk,black_hole.gm_eq[istart:lr,5]);
=#
r0=black_hole.r[56]
α=1.0
sptt = Spline1D(black_hole.r[istart:lr],black_hole.gm_eq[istart:lr,1]); # Needs to discard infinity
sprr = Spline1D(black_hole.r[istart:lr],black_hole.gm_eq[istart:lr,2]);
sptφ = Spline1D(black_hole.r[istart:lr],black_hole.gm_eq[istart:lr,4]);
spφφ = Spline1D(black_hole.r[istart:lr],black_hole.gm_eq[istart:lr,5]);
gm_spline = [sptt,sprr,sptφ,spφφ];
 gtt=evaluate(gm_spline[1],r0);
 gtp=evaluate(gm_spline[3],r0);
  E0=sqrt(-gtt)
  L0=gtp/sqrt(-gtt)
  params=[E0,L0];
   p = (E0, L0, sptt,sprr,sptφ,spφφ,α);
   u₀=[0.0,0.0,r0,0.0]
   prob=ODEProblem(orbit!,u₀,(0.0,10000.0),p)
   @time sol = solve(prob, RK4(),saveat=.1, adaptive=false,dt=0.01);
plot(sol[3,1:istop].*cos.(sol[2,1:istop]),sol[3,1:istop].*sin.(sol[2,1:istop]),xlims=(-1.2,1.2),ylims=(-1.2,1.2),label=false,xlabel=L"x", ylabel=L"y",ratio=1)

gm_spline = [sptt,sprr,sptφ,spφφ];
A = zeros(length(sol.t),2);
sola=zeros(length(sol.t),3);
for i=1:3
    sola[:,i] = sol[i,:]
end
waveform!(A,sola,params,gm_spline; Θ₀=π/4);
tc=time_convertion(black_hole.ADM_M)

plot(sol.t,A[:,1],xlims=(0,500),xlabel=L"t",label=L"A_+")
plot!(sol.t,A[:,2],label=L"A_\times")
=#
function waveform!(A,sol,params,gm_spline; Θ₀=0.0,  Φ₀=0.0)
    E, L = params
    cost = cos(Θ₀)
    Ns = length(sol[:,1])
    for i=1:Ns
        t,p,r,r_τ = sol[i,:]
        δ = p-Φ₀
        ct = cos(δ)
        st = sin(δ)
        c2t = cos(2δ)
        s2t = sin(2δ)

        gtt = evaluate(gm_spline[1],r); gtt_r = derivative(gm_spline[1],r);
        grr = evaluate(gm_spline[2],r); grr_r = derivative(gm_spline[2],r);
        gtp = evaluate(gm_spline[3],r); gtp_r = derivative(gm_spline[3],r);
        gpp = evaluate(gm_spline[4],r); gpp_r = derivative(gm_spline[4],r);

        Δ = gtp^2 - gtt*gpp
        t_τ = (L*gtp + E*gpp)/Δ

        r_t = r_τ/t_τ

        p_t = -(L*gtt + E*gtp)/(L*gtp + E*gpp)

        p_tt = ( - ((L*gtp + E*gpp)*(L*gtt_r + E*gtp_r) - (L*gtp_r + E*gpp_r)*(L*gtt +
        E*gtp)))/(L*gtp + E*gpp)^2*r_t
#=
        r_tt = ((((gtp*gtt_r - 2*gtp_r*gtt)*gtp + gpp_r*gtt^2)*grr - Δ*grr_r*gtt)*L^2 + (
        Δ - E^2*gpp)*Δ*grr_r + (gpp^2*gtt_r - 2*gpp*gtp*gtp_r + gpp_r*gtp^2)*E
        ^2*grr + 2*((gpp_r*gtt - 2*gtp*gtp_r + gpp*gtt_r)*grr*gtp + (grr*gtp_r - grr_r*
        gtp)*Δ)*L*E)/(2*Δ*grr^2*(L*gtp + E*gpp))
=#

        r_tt=r_t*(((gpp_r*gtt - 2*gtp*gtp_r + gpp*gtt_r + E^2*gpp_r + (L*gtt_r + 2*E*gtp_r)*
        L)*(gpp*gtt - gtp^2)*grr*r - (gpp*gtt - gtp^2 + E^2*gpp + (L*gtt + 2*E*
        gtp)*L)*((gpp_r*gtt - 2*gtp*gtp_r + gpp*gtt_r)*grr*r + (gpp*gtt - gtp^2)*(grr
        + grr_r*r)) + 2*(gpp*gtt - gtp^2 + E^2*gpp + (L*gtt + 2*E*gtp)*L)*(gpp_r*
        gtt - 2*gtp*gtp_r + gpp*gtt_r)*grr*r)*(L*gtp + E*gpp) - 2*(gpp*gtt - gtp^2 +
        E^2*gpp + (L*gtt + 2*E*gtp)*L)*(L*gtp_r + E*gpp_r)*(gpp*gtt - gtp^2)*grr
        *r)/(2*sqrt(( - (L*gtt + 2*E*gtp)*L - E^2*gpp - gpp*gtt + gtp^2)/((gpp*gtt
         - gtp^2)*grr*r))*(L*gtp + E*gpp)^2*(gpp*gtt - gtp^2)*grr^2*r^2)

        A[i,1] = (2*(((r*r_tt + r_t^2)*(st + 1)*(st - 1) - (2*st^2 - 1)*p_t^2*r^2 + (4*p_t*
        r_t + p_tt*r)*ct*r*st)*cost^2 + (r*r_tt + r_t^2)*st^2 - (2*st^2 - 1)*p_t^2*
        r^2 + (4*p_t*r_t + p_tt*r)*ct*r*st))


        A[i,2] = (2*(2*(r*r_tt + r_t^2 - 2*p_t^2*r^2)*ct*st - (4*p_t*r_t + p_tt*r)*(2*st^2 -
        1)*r)*cost)

    end
end


function waveform_sol!(A,sol,params; Θ₀=0.0,  Φ₀=0.0)
    E, L = params
    cost = cos(Θ₀)
    sp = Spline1D(sol[:,1],sol[:,2]); sr = Spline1D(sol[:,1],sol[:,3])
    p_ta = derivative(sp,sol[:,1]);
    p_tta = derivative(sp,sol[:,1];nu=2);
    r_ta = derivative(sr,sol[:,1]);
    r_tta = derivative(sr,sol[:,1];nu=2);
    Ns = length(sol[:,1])
    for i=1:Ns
        t,p,r = sol[i,:]
        δ = p-Φ₀
        ct = cos(δ)
        st = sin(δ)
        c2t = cos(2δ)
        s2t = sin(2δ)

        p_t=p_ta[i]; p_tt=p_tta[i];
        r_t=r_ta[i]; r_tt=r_tta[i];


        A[i,1] = (2*(((r*r_tt + r_t^2)*(st + 1)*(st - 1) - (2*st^2 - 1)*p_t^2*r^2 + (4*p_t*
        r_t + p_tt*r)*ct*r*st)*cost^2 + (r*r_tt + r_t^2)*st^2 - (2*st^2 - 1)*p_t^2*
        r^2 + (4*p_t*r_t + p_tt*r)*ct*r*st))

        A[i,2] = (2*(2*(r*r_tt + r_t^2 - 2*p_t^2*r^2)*ct*st - (4*p_t*r_t + p_tt*r)*(2*st^2 -
        1)*r)*cost)

    end
end


#=

io=open("semi_orbit_M=1.1976962626570729.dat","w")
rk=r_to_rk.(rₕ,black_hole.ADM_M,black_hole.ADM_J/black_hole.ADM_M,sol[3,:]);
println(io,"# τ,\t t,\t φ,\t r,\t rk,\t ṙ")
for i=1:length(sol.t)
    println(io,sol.t[i],"\t",sol[1,i],"\t",sol[2,i],"\t",sol[3,i],"\t",rk[i],"\t",sol[4,i])
end
close(io)

io=open("ponty_petal_M=1.1976962626570729.dat","w")
rk=r_to_rk.(rₕ,black_hole.ADM_M,black_hole.ADM_J/black_hole.ADM_M,sol[3,:]);
println(io,"# τ,\t t,\t φ,\t r,\t rk,\t ṙ")
for i=1:length(sol.t)
    println(io,sol.t[i],"\t",sol[1,i],"\t",sol[2,i],"\t",sol[3,i],"\t",rk[i],"\t",sol[4,i])
end
close(io)

istop=1743
semi="../Equatorial_Geodesic/semi_orbit_M=1.1976962626570729.dat"
sols=readdlm(semi,Float64,skipstart=1)
#sol=zeros(length(sols[:,1]),4);
#sol[:,1:3]=sols[:,2:4];sol[:,4]=sols[:,6];
sol=zeros(length(sols[:,1]),3);
sol[:,1:3]=sols[:,2:4];
As=zeros(length(sols[:,1]),2);
r0=sol[1,3];gtt=evaluate(gm_spline[1],r0);gtp=evaluate(gm_spline[3],r0);
E0=sqrt(-gtt);L0=gtp/sqrt(-gtt);

#waveform!(As,sol,[E0,L0],gm_spline);
waveform_sol!(As,sol,[E0,L0]);
io=open("../Equatorial_Geodesic/semi_orbit_M=1.1976962626570729_WF.dat","w")
#rk=r_to_rk.(rₕ,black_hole.ADM_M,black_hole.ADM_J/black_hole.ADM_M,sol[3,:]);
println(io,"# t,\t A₊,\t Aₓ")
for i=1:length(sol[:,1])
    println(io,sol[i,1],"\t",As[i,1],"\t",As[i,2])
end
close(io)

istop=1967
is=2000
pointy="../Equatorial_Geodesic/pointy_petal_M=1.1976962626570729.dat"
solp=readdlm(pointy,Float64,skipstart=1)
#sol=zeros(length(solp[:,1]),4);
#sol[:,1:3]=solp[:,2:4];sol[:,4]=solp[:,6];
sol=zeros(length(solp[:,1]),3);
sol[:,1:3]=solp[:,2:4];
Ap=zeros(length(solp[:,1]),2);
r0=sol[1,3];gtt=evaluate(gm_spline[1],r0);gtp=evaluate(gm_spline[3],r0);
E0=sqrt(-gtt);L0=gtp/sqrt(-gtt);
#waveform!(Ap,sol,[E0,L0],gm_spline);
waveform_sol!(Ap,sol,[E0,L0]);
io=open("../Equatorial_Geodesic/pointy_petal_M=1.1976962626570729_WF.dat","w")
#rk=r_to_rk.(rₕ,black_hole.ADM_M,black_hole.ADM_J/black_hole.ADM_M,sol[3,:]);
println(io,"# t,\t A₊,\t Aₓ")
for i=1:length(sol[:,1])
    println(io,sol[i,1],"\t",Ap[i,1],"\t",Ap[i,2])
end
close(io)
=#
