function edge_exp, x, p0
;because edge-on gals need that wacky bessel function
scalex=abs(x-p0[2])/p0(1)
result=p0(0)*scalex*beselk(scalex,1)
return, result
end
