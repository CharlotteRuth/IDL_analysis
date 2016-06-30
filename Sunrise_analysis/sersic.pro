function sersic, x, p
;p is of form, [bulge csb, bulge Re, bulge n]

;x=abs(x)
;p=abs(p)
N=1

k = 1.9992* p[2] - 0.3271

bulge=p(0)*exp(-k*((x/p(1))^(1/p(2))-1))	;sersic function
result=bulge

return, result
end
