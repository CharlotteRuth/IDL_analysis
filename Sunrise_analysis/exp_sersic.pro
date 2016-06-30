function exp_sersic, x,p
;p=[Disk csb, R_s, bulge csb, R_e, n)


p=abs(p)

disk=abs(p(0))*exp(-x/p(1))	;exponential fit

k = 1.9992* p[4] - 0.3271

bulge=abs(p(2))*exp(-k*((x/p(3))^(1/p(4))-1))     ;sersic fit

result=disk+bulge


return, result
end
