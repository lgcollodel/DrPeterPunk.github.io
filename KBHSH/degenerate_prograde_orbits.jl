using Plots, Glob, DelimitedFiles, Dierckx;

function r(x,rₕ)
    sqrt((x - 1.0)^2*rₕ^2 + x^2)/(1.0-x)
end

function get_metric(isol,rₕ)
    f_read = readdlm(isol,Float64);
    xx = f_read[:,1];
    nt = length(xx);
    for i=1:nt;
        if xx[i]==1.0
            global nx = i
            break
        end
    end
    ny = Int(nt/nx);
    ne = (ny-1)*nx+1;
    x  = xx[ne:nt]
    xp=x[1:end-1];
    y  = f_read[ne:nt,2];
    f0 = f_read[ne:nt,3];
    f1 = f_read[ne:nt,4];
    f2 = f_read[ne:nt,5];
    ω  = -f_read[ne:nt,6];

    gd   = zeros(nx-1,3);
    sqdetg = zeros(nx-1);

    for i=1:nx-1
        rc = r(x[i],rₕ);
        e2f2 = exp(2*f2[i]);
        gd[i,1] = -(1.0-rₕ/rc)*exp(2*f0[i])+e2f2*ω[i]^2;
        gd[i,2] = -e2f2*ω[i]*rc;
        gd[i,3] =  e2f2*rc^2;
        sqdetg[i]  = exp(f0[i]+f2[i]+2*f1[i])*rc^2;
    end
    return xp, gd, sqdetg
end




function get_isco_data(xp,gm,rₕ,n)
    sp1 = Spline1D(xp,gm[:,1]); sp2 = Spline1D(xp,gm[:,2]); sp3 = Spline1D(xp,gm[:,3]);
    ϵ=0.0
    Ω=zeros(n,2)
    en   = zeros(n,2);
    an   = zeros(n,2);
    v_rr = zeros(n,2);
    xi=zeros(n);
    global Lvp_first = false
    global Lvm_first = false
    global Ldiscont  = false
    global L2stable  = false
    for i=2:n
        global Lvp  = false
        global Lvm  = false
        xi[i]=(i-1)/(n+10)
        x0 = xi[i]
        x0 > 0.9 ? break : nothing
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
        if Δ>=0
            Ω[i,1] = (-gtp_r+sqrt(Δ))/gpp_r;
            Ω[i,2] = (-gtp_r-sqrt(Δ))/gpp_r;
            Γ1      = -gtt-2*Ω[i,1]*gtp-Ω[i,1]^2*gpp;
            Γ2      = -gtt-2*Ω[i,2]*gtp-Ω[i,2]^2*gpp;
            if Γ1>0
                en[i,1] = -(gtt+gtp*Ω[i,1])/sqrt(Γ1);
                an[i,1] = (gtp+gpp*Ω[i,1])/sqrt(Γ1);
                v_rr[i,1] = en[i,1]^2*gpp_rr+2*en[i,1]*an[i,1]*gtp_rr+an[i,1]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                if en[i,1] >= 0.0
                    v_rr[i,1] < ϵ ? (global Lvp = true; global Lvp_first = true) : (Lvp_first ? (global Ldiscont = true ) : nothing);
                elseif Lvp_first
                    global Ldiscont=true
                    #println("energy super luminal")
                end
            else
                en[i,1]=NaN; an[i,1]=NaN; v_rr[i,1]=NaN;
                Lvp_first ? (global Ldiscont=true) : nothing
            end
            if Γ2>0
                en[i,2] = -(gtt+gtp*Ω[i,2])/sqrt(Γ2);
                an[i,2] = (gtp+gpp*Ω[i,2])/sqrt(Γ2);
                v_rr[i,2] = en[i,2]^2*gpp_rr+2*en[i,2]*an[i,2]*gtp_rr+an[i,2]^2*gtt_rr-2*gtp_r^2-2*gtp*gtp_rr+gtt_rr*gpp+2*gtt_r*gpp_r+gtt*gpp_rr;
                if (en[i,2] >= 0.0 && v_rr[i,2] <= ϵ && Ω[i,2] >= 0.0)
                    Lvp ? (global L2stable = true) : (global Ldiscont = true)
                end
#=
                if en[i,2] >= 0.0
                    v_rr[i,2] <= ϵ ? (global Lvm = true; ) : nothing;
                end
                (Ω[i,2] > 0.0 && Lvm && Lvp) ? (global L2stable = true) : nothing
=#
            else
                en[i,2]=NaN; an[i,2]=NaN; v_rr[i,2]=NaN;
            end
        else
            Ω[i,1]=NaN; Ω[i,2]=NaN; en[i,1]=NaN; en[i,2]=NaN; an[i,1]=NaN; an[i,2]=NaN; v_rr[i,1]=NaN; v_rr[i,2]=NaN;
            Lvp_first ? (global Ldiscont = true) : nothing
        end
        #println(Ldiscont)
        #Ldiscont ? break : nothing

    end
#    Ldiscont ? id = 1. : id = 0.;
#    L2stable ? i2 = 1. : i2 = 0.;
#    orb_data = [id,i2]
    return  Ldiscont, L2stable
end


cd("/media/Warehouse/Work/SKBH/");
files = glob("n2gc*");

iod = open("discontinuous_prograde_stable_orbits.dat", "a");
io2 = open("degenerate_prograde_stable_orbits.dat", "a");
for ifile in files
#for ifile in files[2:11]
   frh_read = readdlm(ifile, Float64);
#   frh_read = readdlm("om980000.temp",Float64);
    rₕ=0.01; rh="0.01";
    for line=1:size(frh_read)[1]
        nrₕ   = frh_read[line,1];
        if nrₕ!=rₕ
            rₕ = nrₕ;
            length(string(rₕ)) == 4 ? rh = string(rₕ) : rh = string(rₕ,"0");
            Lnewr = true
        else
            Lnewr = false
        end
        ωₛ   = frh_read[line,3]; lom = length(string(ωₛ));
#        om = string(ωₛ)
        lom <= 8 ? om = string(ωₛ,"0"^(8-lom)) : om = string(ωₛ)[1:8];
#        ωₛ=0.980; om=string(ωₛ);
        mass = frh_read[line,4];
        angm = frh_read[line,5];
        surg = frh_read[line,8];
        sol_file = glob("rh=$rh/*$om*")
#        sol_file = glob("om=0.980/*$rh*")#        om = string(ωₛ)

#        println(rh,"\t",om,"\t",sol_file)
        for isol in sol_file
            println(rh,"\t",om,"\t",mass,"\t",isol)
            xx, gd = get_metric(isol,rₕ);

            Ld, L2 = get_isco_data(xx,gd,rₕ,1000);
            orb_data=[rₕ, ωₛ, mass, angm];
            Lnewr ? println(iod,"\n") : nothing
            Lnewr ? println(io2,"\n") : nothing
            Ld ? writedlm(iod,[orb_data]) : nothing
            L2 ? writedlm(io2,[orb_data]) : nothing
#=
            orb_data = get_isco_data(xx,gd,rₕ,1000);
            if  1.0 in orb_data
                Lnewr ? println(io,"\n") : nothing
                pushfirst!(orb_data,rₕ, ωₛ, mass, angm);
                writedlm(io,[orb_data]);
            end
=#
        end
    end
end
close(iod)
close(io2)
