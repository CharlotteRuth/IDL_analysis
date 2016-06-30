function sech_one, x, p0
;fit a sech^2 vertical profile
result=p0(0)*(cosh((x-p0[2])/abs(p0(1))))^(-2.)
return, result
end
