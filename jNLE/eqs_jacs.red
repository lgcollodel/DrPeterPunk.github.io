
eq1 :=  - (exp(f)*(((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2
*sig^4*xlam - sig^4)*dxdr^3*f_x*h_x^2*r + 2*(h + sig)^2*(h - sig)^2*h^2*
xlam)*k1*r - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig
^4*xlam - sig^4)*dxdr_x*r - 2*sig^4)*dxdr*f_x - ((f_x^2 + f_xx)*(h^6*k1*r^
2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4) + 2*(h
+ sig)^2*(h - sig)^2*h^2*h_x^2*k1^2*r^2*xlam)*dxdr^2*r) + ((f_xx + 4*h_x
^2*k1)*dxdr^2*r^2 - 4 + (dxdr_x*r - 3)*dxdr*f_x*r)*h^2*k1*om^2*r*sig^4)

eq2 := exp(f)*(((h^6*h_xx*k1*r^2*xlam - 3*h^5*h_x^2*k1*r^2*xlam - 2*h^4*h_xx
*k1*r^2*sig^2*xlam + 4*h^3*h_x^2*k1*r^2*sig^2*xlam + h^2*h_xx*k1*r^2*sig
^4*xlam - h*h_x^2*k1*r^2*sig^4*xlam - h_xx*sig^4 + (h^6*k1*r^2*xlam - 2*h
^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*f_x*h_x)*dxdr^2 -
((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam -
sig^4)*dxdr^3*h_x^3*k1*r - (3*h^2 - sig^2)*(h + sig)*(h - sig)*h*xlam))*r +
 (2*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam -
 sig^4)*h_x + (3*h^2 - sig^2)*(h + sig)*(h - sig)*f_x*h*r^2*xlam + (h^6*k1*
r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*
dxdr_x*h_x*r)*dxdr) - ((f_x + h*h_x*k1 + dxdr_x*h*h_x*k1*r)*dxdr*r + (h*h_xx -
h_x^2)*dxdr^2*k1*r^2 + 1)*h*om^2*r*sig^4

ju(1,1) :=  - exp(f)*(((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r
^2*sig^4*xlam - sig^4)*dxdr^3*f_x*h_x^2*r + 2*(h + sig)^2*(h - sig)^2*h^
2*xlam)*k1*r - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*
sig^4*xlam - sig^4)*dxdr_x*r - 2*sig^4)*dxdr*f_x - ((f_x^2 + f_xx)*(h^6*k1*
r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4) + 2*
(h + sig)^2*(h - sig)^2*h^2*h_x^2*k1^2*r^2*xlam)*dxdr^2*r)

jux(1,1) := (exp(f)*((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
2*sig^4*xlam - sig^4)*dxdr_x*r - 2*sig^4 - (dxdr*h_x^2*k1*r - 2*f_x)*(h^6*
k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*
dxdr*r) - (dxdr_x*r - 3)*h^2*k1*om^2*r^2*sig^4)*dxdr

juxx(1,1) := (exp(f)*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
2*sig^4*xlam - sig^4) - h^2*k1*om^2*r^2*sig^4)*dxdr^2*r

ju(1,2) :=  - 2*(exp(f)*((dxdr^2*h_x^2*k1*r - dxdr_x)*dxdr*f_x*r^2 + 2 - (f_x^
2 + f_xx + 2*h_x^2*k1)*dxdr^2*r^2)*(3*h^2 - sig^2)*(h + sig)*(h - sig)*xlam
 + ((f_xx + 4*h_x^2*k1)*dxdr^2*r^2 - 4 + (dxdr_x*r - 3)*dxdr*f_x*r)*om^2*sig
^4)*h*k1*r

jux(1,2) :=  - 2*(exp(f)*((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*
k1*r^2*sig^4*xlam - sig^4)*dxdr*f_x - 2*(h + sig)^2*(h - sig)^2*h^2*k1*r*
xlam) + 4*h^2*k1*om^2*r*sig^4)*dxdr^2*h_x*k1*r^2

juxx(1,2) := 0

ju(2,1) := exp(f)*(((h^6*h_xx*k1*r^2*xlam - 3*h^5*h_x^2*k1*r^2*xlam - 2*h^4*
h_xx*k1*r^2*sig^2*xlam + 4*h^3*h_x^2*k1*r^2*sig^2*xlam + h^2*h_xx*k1*r^2
*sig^4*xlam - h*h_x^2*k1*r^2*sig^4*xlam - h_xx*sig^4 + (h^6*k1*r^2*xlam -
 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*f_x*h_x)*dxdr^
2 - ((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
- sig^4)*dxdr^3*h_x^3*k1*r - (3*h^2 - sig^2)*(h + sig)*(h - sig)*h*xlam))*r
 + (2*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam
 - sig^4)*h_x + (3*h^2 - sig^2)*(h + sig)*(h - sig)*f_x*h*r^2*xlam + (h^6*
k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*
dxdr_x*h_x*r)*dxdr)

jux(2,1) := (exp(f)*((h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
2*sig^4*xlam - sig^4)*dxdr*h_x + (3*h^2 - sig^2)*(h + sig)*(h - sig)*h*r*
xlam) - h*om^2*r*sig^4)*dxdr*r

juxx(2,1) := 0

ju(2,2) := (exp(f)*((6*h^5*h_xx - 15*h^4*h_x^2 - 8*h^3*h_xx*sig^2 + 12*h^2*
h_x^2*sig^2 + 2*h*h_xx*sig^4 - h_x^2*sig^4 + 2*(3*h^2 - sig^2)*(h + sig)*
(h - sig)*f_x*h*h_x)*dxdr^2*k1*r^2 + 15*h^4 - 12*h^2*sig^2 + sig^4 - 2*(3*
h^2 - sig^2)*(h + sig)*(h - sig)*dxdr^3*h*h_x^3*k1^2*r^3 + (2*(dxdr_x*r +
2)*(3*h^2 - sig^2)*(h + sig)*(h - sig)*h*h_x*k1 + (15*h^4 - 12*h^2*sig^2 +
sig^4)*f_x)*dxdr*r)*xlam - ((f_x + 2*h*h_x*k1 + 2*dxdr_x*h*h_x*k1*r)*dxdr*r + (
2*h*h_xx - h_x^2)*dxdr^2*k1*r^2 + 1)*om^2*sig^4)*r

jux(2,2) := (exp(f)*((dxdr_x*r + 2 - 3*dxdr^2*h_x^2*k1*r^2)*(h^6*k1*r^2*xlam
- 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4) + ((h^6*k1*r
^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^2*sig^4*xlam - sig^4)*f_x -
 2*(3*h^2 - sig^2)*(h + sig)*(h - sig)*h*h_x*k1*r^2*xlam)*dxdr*r) - ((dxdr_x*
r + 1)*h - 2*dxdr*h_x*r)*h*k1*om^2*r^2*sig^4)*dxdr

juxx(2,2) := (exp(f)*(h^6*k1*r^2*xlam - 2*h^4*k1*r^2*sig^2*xlam + h^2*k1*r^
2*sig^4*xlam - sig^4) - h^2*k1*om^2*r^2*sig^4)*dxdr^2*r
