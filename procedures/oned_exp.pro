function oned_exp, x, p0

result=p0(0)*exp(-1.*abs(x-p0[2])/p0(1))

return, result
end


