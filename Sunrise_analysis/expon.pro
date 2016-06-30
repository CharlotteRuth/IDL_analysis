function expon, x, p
;p is of form, [disk csb, disk Rl]
x=abs(x)
p=abs(p)
N=1


disk=p(0)*exp(-x/p(1))	;exponential fit
result=disk

return, result
end
