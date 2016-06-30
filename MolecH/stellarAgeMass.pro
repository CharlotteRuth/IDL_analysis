; .r /astro/users/christensen/code/MolecH/envelope_LW.pro

pro stellarAgeMass
loadct,39
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/'
file = 'h516.cosmo25cmb.1536g1MBWK.00512'
params = tipsyunits('/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g1MBWK/h516.cosmo25cmb.1536g1MBWK.param')

rtipsy,dir+file,h,g,d,s
readcol,dir+file+'.massform',ms2
ms = ms2[N_ELEMENTS(ms2) - h.nstar:N_ELEMENTS(ms2) - 1]
window,0
plot,(MAX(s.tform) - s.tform)*params.timeunit,s.mass/ms,psym = 3,/xlog,/ylog,yrange = [0.7,1.0],ystyle = 1,xrange = [1e6,1e10]

age = (MAX(s.tform) - s.tform)*params.timeunit
ind = where(age gt 1e6 AND age lt 5e9) 
age = age[ind]
massfrac = s[ind].mass/ms[ind]
oplot,age,massfrac,psym = 3,color = 240
stop

logage = alog10(age*31556926)
logmassfrac = alog10(massfrac)
fit = POLY_FIT(logmassfrac,logage,4)
logmassfrac_fit = (findgen(100)*(0 + 0.14)/100.0 -0.14)
logage_fit =  fit[0] + $
              fit[1]*logmassfrac_fit + $
              fit[2]*logmassfrac_fit*logmassfrac_fit + $ 
              fit[3]*logmassfrac_fit*logmassfrac_fit*logmassfrac_fit + $
              fit[4]*logmassfrac_fit*logmassfrac_fit*logmassfrac_fit*logmassfrac_fit
window,1
plot,10^logmassfrac_fit,10^logage_fit,/ylog,/xlog,xrange=[0.7,1],xstyle = 1
oplot,massfrac,age,psym = 3,color = 240

print,fit
stop

end
