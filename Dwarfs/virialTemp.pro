;virialTemp,'/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/1E6R/10M_k/o10M_1.00300'
PRO virialTemp,filename
; Calculates the virial temperature, as in BT
C = 6.8d5
mu = 0.61
scaleM = 1e11
scaleR = 10

lengthunit = 3.0857d21 ;system length unit in cm (=1kpc)
kpc_per_syslength = 1.0
msol_per_sysmass = 2.362d5 ;system mass unit in solar masses
nbins = 500

rtipsy,filename,h,g,d,s
allmass = msol_per_sysmass*[g.mass,d.mass,s.mass]
totalmass = TOTAL(allmass)
allr = kpc_per_syslength*[SQRT(g.x*g.x + g.y*g.y + g.z*g.z),SQRT(d.x*d.x + d.y*d.y + d.z*d.z),SQRT(s.x*s.x + s.y*s.y + s.z*s.z)]
sort_r = SORT(allr)
allmass = allmass[sort_r]
allr = allr[sort_r]
histr = histogram(allr,nbins = nbins,reverse_indices = R, locations = xbins,max = 100)
histMass = fltarr(nbins)
histMassTotal = fltarr(nbins)
histMass[0] = TOTAL(allmass[R[R[0] : R[0+1]-1]])
histMassTotal[0] = TOTAL(allmass[R[R[0] : R[0+1]-1]])
FOR i = 1,nbins - 1 DO BEGIN
    IF (R[i] NE R[i+1]) THEN histMass[i] = TOTAL(allmass[R[R[I] : R[i+1]-1]]) + histMass[i - 1] ELSE histMass[i] = 0
    IF (R[i] NE R[i+1]) THEN histMassTotal[i] = TOTAL(allmass[R[R[I] : R[i+1]-1]]) + histMassTotal[i - 1] ELSE histMassTotal[i] = histMassTotal[i - 1]
ENDFOR
set_plot,'x'
plot,xbins,histMassTotal/totalmass,/xlog,xrange = [1,100]
oplot,[1,800],[0.5,0.5],linestyle = 2
temp = MIN(ABS(histMassTotal/totalmass - 0.5),ind)
print,'r_h: ',xbins[ind]
oplot,[xbins[ind],xbins[ind]],[0,1],linestyle = 2
T_vir = C*mu*totalmass/scaleM*scaleR/xbins[ind]
print,T_vir
stop
END
