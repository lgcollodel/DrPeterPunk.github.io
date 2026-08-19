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
#        sqdetg[i]  = exp(f0[i]+f2[i]+2*f1[i])*rc^2;
        sqdetg[i]  = exp(f0[i]+f2[i]+f1[i])*rc;

    end
    return xp, gd, sqdetg
end
function get_metric_bs(isol)
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
    f = f_read[ne:nt,3];
    l = f_read[ne:nt,4];
    g = f_read[ne:nt,5];
    ω  = -f_read[ne:nt,6];

    gd   = zeros(nx-1,3);
    sqdetg = zeros(nx-1);

    for i=1:nx-1
        rc = r(x[i],0.0);
        gd[i,1] = -f[i]+l[i]/f[i]*ω[i]^2;
        gd[i,2] = -l[i]/f[i]*ω[i]*rc;
        gd[i,3] =  l[i]/f[i]*rc^2;
        sqdetg[i]  = (sqrt(l[i])*g[i]*l[i]*rc^2)/f[i]
    end
    return xp, gd, sqdetg
end
