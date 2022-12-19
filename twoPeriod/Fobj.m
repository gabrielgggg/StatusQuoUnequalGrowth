function f = Fobj(x)

global k beta p iui2 iuj2

f = -(log(x(1))+k*log(x(3))+...
    beta*(p*iui2(x(2))+(1-p)*iuj2(x(1))));