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

function lightring_radius(xp,gm,rₕ,n)
    sp1 = Spline1D(xp,gm[:,1]); sp2 = Spline1D(xp,gm[:,2]); sp3 = Spline1D(xp,gm[:,3]);
    v    = zeros(n,2)
    v_r  = zeros(n,2);
    xi=zeros(n);
    global lr₊ = false
    global lr₋ = false
    for i=2:n
        global Lvp  = false
        global Lvm  = false
        xi[i]=(i-1)/(n+10)
        x0 = xi[i]
        x0 > 0.9 ? break : nothing
        dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
        gtt     = evaluate(sp1,x0);
        gtt_r   = derivative(sp1,x0)*dx_dr;
        gtp     = evaluate(sp2,x0);
        gtp_r   = derivative(sp2,x0)*dx_dr;
        gpp     = evaluate(sp3,x0);
        gpp_r   = derivative(sp3,x0)*dx_dr;
        Δ       = gtp^2-gtt*gpp;
        if Δ>=0
            vb       = sqrt(Δ)
            v[i,1]   = ( -gtp + vb ) / gpp
            v[i,2]   = ( -gtp - vb ) / gpp
            v_r[i,1] = ( -gtp_r*gpp+gtp*gpp_r )/gpp^2 + ( 2gtp*gtp_r*gpp - gtt_r*gpp^2 + gtt*gpp_r*gpp - 2*gpp_r*gtp^2)/(2*gpp^2*vb)
            v_r[i,2] = ( -gtp_r*gpp+gtp*gpp_r )/gpp^2 - ( 2gtp*gtp_r*gpp - gtt_r*gpp^2 + gtt*gpp_r*gpp - 2*gpp_r*gtp^2)/(2*gpp^2*vb)
            (!lr₊ && v_r[i,1]*v_r[i-1,1] < 0.0) ? (global lr₊ = true; global i_lr₊=i;) : nothing;
            (!lr₋ && v_r[i,2]*v_r[i-1,2] < 0.0) ? (global lr₋ = true; global i_lr₋=i;) : nothing;

        else
            v[i,1]=NaN; v[i,2]=NaN; v_r[i,1]=NaN; v_r[i,2]=NaN; a
        end

        (lr₊ && lr₋) ? break : nothing

    end
    lightring = [ xi[i_lr₊], v[i_lr₊,1], xi[i_lr₋], v[i_lr₋,2]]
    return lightring
end





cd("/media/Warehouse/Work/SKBH/");
files = glob("n2gc*");
io = open("lightrings.dat", "a");
#for ifile in files[22:30]
for ifile in files[2:2]
   frh_read = readdlm(ifile, Float64);
#   frh_read = readdlm("om980000.temp",Float64);
    rₕ=0.01; rh="0.01";
    for line=1:size(frh_read)[1]
        nrₕ   = frh_read[line,1];
        if nrₕ!=rₕ
            rₕ = nrₕ;
            length(string(rₕ)) == 4 ? rh = string(rₕ) : rh = string(rₕ,"0");
            println(io,"\n")
        end
        ωₛ   = frh_read[line,3]; lom = length(string(ωₛ));
#        om = string(ωₛ)
        lom <= 8 ? om = string(ωₛ,"0"^(8-lom)) : om = string(ωₛ)[1:8];
#        ωₛ=0.980; om=string(ωₛ);
        mass = frh_read[line,4];
        angm = frh_read[line,5];
        surg = frh_read[line,8];
        sol_file = glob("rh=$rh/*$om*")
#        sol_file = glob("om=0.980/*$rh*")
#        println(rh,"\t",om,"\t",sol_file)
        for isol in sol_file
            println(rh,"\t",om,"\t",mass,"\t",isol)
            xx, gd = get_metric(isol,rₕ);
            lightring = lightring_radius(xx, gd, rₕ, 1000)
            pushfirst!(lightring, rh, om,  mass, angm);
            writedlm(io,[lightring]);
        end
    end
end
close(io)
