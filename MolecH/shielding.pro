function S_H2,yH2,h
omega_H2 = 0.2
x = yH2*h/5d14
return, (1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))
end

function S_d,yHI,yH2,h,z
sigma_d = 2d-21
zsol = 0.025
return, exp(-1.0*sigma_d*(yHI*h + 2.0*yH2*h))
end
