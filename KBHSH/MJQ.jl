using DelimitedFiles, Dierckx, Glob



function charge_densities(isol,rₕ,ωₛ)
    m=1
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
     θ = zeros(ny)

    gdₕ = zeros(ny,5)
    dmv = zeros(nx,ny)
    djv = zeros(nx,ny)
    dmh = zeros(ny)
    djh = zeros(ny)
      o = zeros(nx,ny)
    ne = (ny-1)*nx+1;
    x  = xx[ne:nt]
    for j=1:ny
        θ[j] = f_read[nx*j,2]
         θ0  = θ[j]
        for i=1:nx
            x0 = x[i]
            r0 =  sqrt((x0 - 1.0)^2*rₕ^2 + x0^2)/(1.0-x0)
            dxdr   = r0*(1.0-x0)^3/x0;
            f0 =  f_read[i+(j-1)*nx,3]
            f1 =  f_read[i+(j-1)*nx,4]
            f2 =  f_read[i+(j-1)*nx,5]
            o[i,j] =  f_read[i+(j-1)*nx,6]
            ω = o[i,j]
             ψ =  f_read[i+(j-1)*nx,7]
            if i==1
                gdₕ[j,1] = f0
                gdₕ[j,2] = f1
                gdₕ[j,3] = f2
                gdₕ[j,4] = ω
                gdₕ[j,5] = ψ
            end

            if (i==1 || i==nx || j==1)
                dmv[i,j]=0.0
                djv[i,j]=0.0
            else
                dmv[i,j] = exp(2*f1-f0+f2)*(-2*(ω*m+ωₛ*r0)*ωₛ + (r0-rₕ)*exp(2*f0))*ψ^2*r0^2*sin(θ0)/((r0-rₕ)*dxdr)
                djv[i,j] = exp(2*f1-f0+f2)*(ω*m+ωₛ*r0)*ψ^2*r0^2*sin(θ0)/((r0-rₕ)*dxdr)
            end
             #=
             for k =1:5
                 gd[i,j,k] = f_read[i+(j-1)*nx,k+2]
             end
             =#
        end
    end

    dmy = zeros(ny)
    djy = zeros(ny)
    for j=1:ny
        θ0 = θ[j]
        spo = Spline1D(x,o[:,j])
        ωxx = derivative(spo,0,nu=2)
        f0  = gdₕ[j,1]
        f1  = gdₕ[j,2]
        f2  = gdₕ[j,3]
         ω  = gdₕ[j,4]
         ψ  = gdₕ[j,5]

         dmh[j] = 0.5*exp(f2-f0)*((ω-ωxx*rₕ^2)*exp(2*f2)*ω*sin(θ0)^2+exp(2*f0))*rₕ*sin(θ0)
         djh[j] = -0.5*exp(f2-f0)*(ω-ωxx*rₕ^2)*exp(2*f2)*rₕ^2*sin(θ0)^3

         dmx = dmv[:,j]
         spm = Spline1D(x,dmx)
         dmy[j] = integrate(spm,0,1)

         djx = djv[:,j]
         spj = Spline1D(x,djx)
         djy[j] = integrate(spj,0,1)

    end

    smv = Spline1D(θ, dmy)
    sjv = Spline1D(θ, djy)
    smh = Spline1D(θ, dmh)
    sjh = Spline1D(θ, djh)

    mv = -integrate(smv, 0, π/2)
    jv =  integrate(sjv, 0, π/2)
    mh =  integrate(smh, 0, π/2)
    jh =  integrate(sjh, 0, π/2)/2

    return x, θ, mv, jv, mh, jh
end
