function pot_free(h,p)
    k1,xlam,om,sig = p
    pot    = xlam/2*h^2
    pot_h  = xlam*h
    pot_hh = xlam
    return pot,pot_h,pot_hh
end
function pot_solitonic(h,p)
    k1,xlam,om,sig = p
    pot    = xlam*h^2*(1.0-h^2/sig^2)^2
    pot_h  = 2*pot/h-4*pot/(1.0-h^2/sig^2)*h/sig^2
    pot_hh = (2*(15*h^4 - 12*h^2*sig^2 + sig^4)*xlam)/sig^4
    return pot,pot_h,pot_hh
end
function ueqs!(Pu,u,ux,uxx,p,x)
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
function ujac_eq!(Ju,Jux,Juxx,u,ux,uxx,p,x)
    k1,xlam,om,sig = p
    f, h = u
    f_x,h_x = ux
    f_xx,h_xx = uxx
    r = x/(1.0-x)
    dxdr = (1.0-x)^2
    dxdr_x = 2*(x-1.0)

    pot, pot_h, pot_hh = potential(h,p)

    Ju[1,1] =  - exp(f)*(((k1*pot*r^2 - 1)*dxdr^3*f_x*h_x^2*r + 2*pot)*k1*r - ((k1
    *pot*r^2 - 1)*dxdr_x*r - 2)*dxdr*f_x - ((f_x^2 + f_xx)*(k1*pot*r^2 - 1) + 2*
    h_x^2*k1^2*pot*r^2)*dxdr^2*r)

    Jux[1,1] = (exp(f)*((k1*pot*r^2 - 1)*dxdr_x*r - 2 - (dxdr*h_x^2*k1*r - 2*f_x)*(
    k1*pot*r^2 - 1)*dxdr*r) - (dxdr_x*r - 3)*h^2*k1*om^2*r^2)*dxdr

    Juxx[1,1] = (exp(f)*(k1*pot*r^2 - 1) - h^2*k1*om^2*r^2)*dxdr^2*r

    Ju[1,2] =  - (2*((f_xx + 4*h_x^2*k1)*dxdr^2*r^2 - 4 + (dxdr_x*r - 3)*dxdr*
    f_x*r)*h*om^2 - exp(f)*((f_xx + 2*h_x^2*k1 + f_x^2)*dxdr^2*r^2 - (dxdr^3*f_x
    *h_x^2*k1*r^3 - dxdr*dxdr_x*f_x*r^2 + 2))*pot_h)*k1*r

    Jux[1,2] =  - 2*(exp(f)*((k1*pot*r^2 - 1)*dxdr*f_x - 2*k1*pot*r) + 4*h^2*k1*om
    ^2*r)*dxdr^2*h_x*k1*r^2

    Juxx[1,2] = 0

    Ju[2,1] = exp(f)*(((2*(f_x*h_x + h_xx)*(k1*pot*r^2 - 1) - h_x^2*k1*pot_h*r^2)*
    dxdr^2 - (2*(k1*pot*r^2 - 1)*dxdr^3*h_x^3*k1*r - pot_h))*r + (2*(dxdr_x*r +
    2)*(k1*pot*r^2 - 1)*h_x + f_x*pot_h*r^2)*dxdr)

    Jux[2,1] = (exp(f)*(2*(k1*pot*r^2 - 1)*dxdr*h_x + pot_h*r) - 2*h*om^2*r)*dxdr*r


    Juxx[2,1] = 0

    Ju[2,2] =  - (exp(f)*((h_x^2*pot_hh - 2*h_xx*pot_h - 2*f_x*h_x*pot_h)*dxdr^2*k1
    *r^2 + 2*dxdr^3*h_x^3*k1^2*pot_h*r^3 - pot_hh - (f_x*pot_hh + 4*h_x*k1*
    pot_h + 2*dxdr_x*h_x*k1*pot_h*r)*dxdr*r) + 2*((f_x + 2*h*h_x*k1 + 2*dxdr_x*h*h_x
    *k1*r)*dxdr*r + (2*h*h_xx - h_x^2)*dxdr^2*k1*r^2 + 1)*om^2)*r

    Jux[2,2] = 2*(exp(f)*((dxdr_x*r + 2 - 3*dxdr^2*h_x^2*k1*r^2)*(k1*pot*r^2 - 1)
     + ((k1*pot*r^2 - 1)*f_x - h_x*k1*pot_h*r^2)*dxdr*r) - ((dxdr_x*r + 1)*h - 2*
    dxdr*h_x*r)*h*k1*om^2*r^2)*dxdr

    Juxx[2,2] = 2*(exp(f)*(k1*pot*r^2 - 1) - h^2*k1*om^2*r^2)*dxdr^2*r
end



function ubc!(c,u,ux,p)
    f, h = u
    f_x,h_x = ux
    k1,xlam,om,sig = p
    # origin
    c[1]=f_x
    c[2]=h_x
    # infinity
    c[3]=f
    c[4]=h
end
function ujac_bc!(dc,dcx,u,ux)
    f, h = u
    f_x,h_x = ux
    # origin
    dcx[1,1]=1.0
    dcx[2,2]=1.0
    # infinity
    dc[3,1]=1.0
    dc[4,2]=1.0
end
