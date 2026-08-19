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
    xi=zeros(n);
    global Lfirst = false
    global Ldiscont  = false
    for i=2:n
        global Lvp  = false
        global Lvm  = false
        xi[i]=(i-1)/(n+10)
        x0 = xi[i]
        x0 > 0.9 ? break : nothing
        dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
        gtt_r   = derivative(sp1,x0)*dx_dr;
        gtp_r   = derivative(sp2,x0)*dx_dr;
        gpp_r   = derivative(sp3,x0)*dx_dr;
        Δ       = gtp_r^2-gtt_r*gpp_r;
        if Δ>=0
            global Lfirst = true
        else
            Lfirst ? (global Ldiscont = true; break;) : nothing
        end
    end
    return  Ldiscont
end


cd("/media/Warehouse/Work/SKBH/");
files = glob("n2gc*");
io = open("discontinuous_omega2.dat", "a");
for ifile in files
#for ifile in files[2:31]
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
            Ldiscont = get_isco_data(xx,gd,rₕ,1000);
            if  Ldiscont
                orb_data = [rₕ, ωₛ, mass, angm];
                writedlm(io,[orb_data]);
            end
        end
    end
end
close(io)
