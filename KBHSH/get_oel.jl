
function get_oel(xp,gm,rₕ,sqdetg,n)
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
    Ωp[1]=NaN; Ωm[1]=NaN;Ω[1]=NaN;v2p[1]=NaN;v2m[1]=NaN; Ωsp[1]=NaN; Ωsm[1]=NaN;
#    v_rr = zeros(n,2);
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
        if Δ>=0
            Ωp[i] = (-gtp_r+sqrt(Δ))/gpp_r;
            Ωm[i] = (-gtp_r-sqrt(Δ))/gpp_r;
#            (!Ldegenerated && Ω[i,2] > 0.0) ? global Ldegenerated=true : nothing;
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
#        (Lvp && Lvm) ? println(i, "\t",xi[i],"\t", Ωp[i], "\t", Ωm[i]) : nothing

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
#    return xi, Ω, Ωp, Ωm, en, an, a₊, a₋, v2p, v2m, Δr, sqdg
#    return xi, Ωp, Ωm, Ωsp, Ωsm
    return xi, Ω,  en, an, sqdg

end
