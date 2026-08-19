function eqs1(Pu,u,ux,uxx,p,x)
    k1,xlam,om,sig = p
    f, h = u
    f_x,h_x = ux
    f_xx,h_xx = uxx
    r = x/(1.0-x)
    dxdr = (1.0-x)^2
    dxdr_x = 2*(x-1.0)


    Pu[1] =  - (exp(f)*(((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r
    ^2*sig^4*xlam - sig^4)*dxdr^3*f_x*h_x^2*r + 2*(h + sig)^2*(h - sig)^2*h^
    2*xlam)*k1*r - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*
    sig^4*xlam - sig^4)*dxdr_x*r - 2*sig^4)*dxdr*f_x - ((f_x^2 + f_xx)*(h^6*k1*
    r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4) + 2*
    (h + sig)^2*(h - sig)^2*h^2*h_x^2*k1^2*r^2*xlam)*dxdr^2*r) + ((f_xx + 4*
    h_x^2*k1)*dxdr^2*r^2 - 4 + (dxdr_x*r - 3)*dxdr*f_x*r)*h^2*k1*om^2*r*sig^4)


    Pu[2] = exp(f)*(((h^6*h_xx*k1*r^2*xlam - 3*h^5*h_x^2*k1*r^2*xlam - 2*h^4*
    h_xx*k1*r^2*sig^2*xlam + 4*h^3*h_x^2*k1*r^2*sig^2*xlam + h^2*h_xx*k1*r^2
    *sig^4*xlam - h*h_x^2*k1*r^2*sig^4*xlam - h_xx*sig^4 + (h^6*k1*r^2*xlam -
     2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*f_x*h_x)*dxdr^
    2 - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
    - sig^4)*dxdr^3*h_x^3*k1*r - (3*h^2 - sig^2)*(h + sig)*(h - sig)*h*xlam))*r
     + (2*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
     - sig^4)*h_x + (3*h^2 - sig^2)*(h + sig)*(h - sig)*f_x*h*r^2*xlam + (h^6*
    k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*
    dxdr_x*h_x*r)*dxdr) - ((f_x + h*h_x*k1 + dxdr_x*h*h_x*k1*r)*dxdr*r + (h*h_xx -
    h_x^2)*dxdr^2*k1*r^2 + 1)*h*om^2*r*sig^4
end

function eqs2(Pu,u,ux,uxx,p,x)
    k1,xlam,om,sig = p
    f, h = u
    f_x,h_x = ux
    f_xx,h_xx = uxx
    r = x/(1.0-x)
    dxdr = (1.0-x)^2
    dxdr_x = 2*(x-1.0)

    pot, pot_h, pot_hh = potential(h,p)

    Pu[1] =  - (exp(f)*(((k1*pot*r^2 - 1)*dxdr^3*f_x*h_x^2*r + 2*pot)*k1*r - ((k1*
    pot*r^2 - 1)*dxdr_x*r - 2)*dxdr*f_x - ((f_x^2 + f_xx)*(k1*pot*r^2 - 1) + 2*
    h_x^2*k1^2*pot*r^2)*dxdr^2*r) + ((f_xx + 4*h_x^2*k1)*dxdr^2*r^2 - 4 + (
    dxdr_x*r - 3)*dxdr*f_x*r)*h^2*k1*om^2*r)

    Pu[2] = exp(f)*(((2*(f_x*h_x + h_xx)*(k1*pot*r^2 - 1) - h_x^2*k1*pot_h*r^2)*dxdr
    ^2 - (2*(k1*pot*r^2 - 1)*dxdr^3*h_x^3*k1*r - pot_h))*r + (2*(dxdr_x*r + 2)*(
    k1*pot*r^2 - 1)*h_x + f_x*pot_h*r^2)*dxdr) - 2*((f_x + h*h_x*k1 + dxdr_x*h*h_x
    *k1*r)*dxdr*r + (h*h_xx - h_x^2)*dxdr^2*k1*r^2 + 1)*h*om^2*r
end

for i=1:dform.Nx
    eqs1(Pu,u0[i,:],vu.ux[i,:],vu.uxx[i,:],p,x[i])
    eqs2(Pv,u0[i,:],vu.ux[i,:],vu.uxx[i,:],p,x[i])
    diff[i,1]=abs(Pu[1]-Pv[1])
    diff[i,2]=abs(Pu[2]-Pv[2])
end

# Normalization is different
for i=1:dform.Nx
           eqs1(Pu,u0[i,:],vu.ux[i,:],vu.uxx[i,:],p,x[i])
           eqs2(Pv,u0[i,:],vu.ux[i,:],vu.uxx[i,:],p,x[i])
           diff[i,1]=abs(Pu[1]-Pv[1]*0.1^4)
           diff[i,2]=abs(Pu[2]-Pv[2]*0.1^4)
end



# Try fine grid
dx=zeros(202);
for i=1:202
    dx[i]=x0[i+1]-x0[i]
end
dx=findmin(dx)[1]

Nx=4109;

u1=zeros(length(x0),2);
for i=1:length(x0)
       u1[i,1]=newvaluate(x,u0[:,1],Px,x0[i])
       u1[i,2]=newvaluate(x,u0[:,2],Px,x0[i])
end
dx=x[2]-x[1]
dxmin=dx/2
0.00024342602314197936
dxmin=dx/2
0.00012171301157098968
xmin=0.9;xmax=0.975
x1=collect(range(0,xmin,length=380))
x2=collect(range(xmin+dxmin,xmax,length=500))
x3=collect(range(xmax+dxmin,1,length=40))
x0=vcat(x1,x2,x3)

#save sol u00=u0;x00=x;
x=x0;u0=u1;

i1=1
xmin=0.97;xmax=1;
for i=1:length(x)
    x[i]>=xmin ? (i1=i-1;break) : nothing
end
Nxn=i1+5*(length(x)-i1)
x0=x[1:i1]
for i=i1:length(x)-1
    push!(x0,x0[end]+(x[i+1]-x[i])/2)
    push!(x0,x[i+1])
end
# Plots

using Glob, Dierckx, DelimitedFiles, Plots
files=glob("h0*sigma=0.100");
Nf=length(files);
oma=[];Ma=[];ha=[];
for i=1:Nf
    file=files[i]
    i0=findfirst("om=",file)[end]+1
    ie=findfirst("_sigma",file)[1]-1
    om1=parse(Float64,file[i0:ie])
    s0 = readdlm(file,Float64);
    x=s0[:,1]; f=s0[:,2]; h=s0[:,3];
    sf=Spline1D(x,f);
    push!(oma,om1)
    push!(Ma,derivative(sf,1)/2)
    push!(ha,h[1])
end
