


function get_light_rings(xp,gm,rₕ,n)
    o   = -gm[:,2]./gm[:,3]
    Vp  = o-sqrt.(o.^2-gm[:,1]./gm[:,3])
    Vm  = o+sqrt.(o.^2-gm[:,1]./gm[:,3])
    spp = Spline1D(xp,Vp);
    spm = Spline1D(xp,Vm);


    lrp = []
    lrm = []
    vp_r = zeros(n)
    vm_r = zeros(n)
    vp_r[1] = NaN
    vm_r[1] = NaN

    for i=2:n
        x0 = (i-1)/(n+10)
        x0 > 0.9 ? break : nothing
        dx_dr   = sqrt((x0 - 1.)^2*rₕ^2 + x0^2)*(x0 - 1.)^2/x0;
        vp = evaluate(spp,x0)
        vm = evaluate(spm,x0)
        vp_r[i] = derivative(spp,x0)*dx_dr
        vm_r[i] = derivative(spm,x0)*dx_dr
        vp_r[i-1]*vp_r[i] <= 0 ? push!(lrp,(x0,vp)) : nothing
        vm_r[i-1]*vp_r[i] <= 0 ? push!(lrm,(x0,vm)) : nothing
    end

    return lrp, lrm

end
