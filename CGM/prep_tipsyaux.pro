;filename = 'h277.cosmo50cmb.3072g14HMbwK.00512.1'
;pfile = '../../h277.cosmo50cmb.3072g14HMbwK.param' 

PRO fix_smooth,filename,filename2 = filename2
IF NOT keyword_set(filename2) THEN filename2 = filename
rtipsy,filename + '.std',h,g,d,s
readarr,filename2 + '.hsmooth',h,smooth,/ascii,part = 'gas'
g.h = smooth
wtipsy,filename + '.std',h,g,d,s,/standard
END

PRO prep_tipsyaux,filename,pfile
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
grav = 6.6725985e-8 ;am^3/g/s
s_per_yr = 3.15569e7
Carbon   = 0.00213
Oxygen   = 0.00541
Silicon  = 0.00067
Iron     = 0.00117

units = tipsyunits(pfile)
rtipsy,filename + '.std',h,g,d,s
dsfr   = fltarr(h.ndark)
ssfr   = fltarr(h.nstar)
gsfr   = fltarr(h.ngas)
dummyf = fltarr(h.ngas)
dummyi = lonarr(h.ngas) + 826393640
readarr,filename + '.H2'        ,h,H2  ,/ascii,part = 'gas'
readarr,filename + '.HI'        ,h,HI  ,/ascii,part = 'gas'
readarr,filename + '.HeI'       ,h,HeI ,/ascii,part = 'gas'
readarr,filename + '.HeII'      ,h,HeII,/ascii,part = 'gas'
readarr,filename + '.OxMassFrac',h,OX  ,/ascii,part = 'gas'
readarr,filename + '.FeMassFrac',h,FE  ,/ascii,part = 'gas'

cstar = 0.1
tempcut = 1e3
denscut = 0.1
deltat = 1e6*s_per_yr ;yr
indsf = where(g.tempg LE tempcut AND g.dens*units.rhounit GE denscut)
tdyn = 1.0/sqrt(4.0*!PI*grav*g.dens*units.rhounit*gm_per_H)
gsfr[indsf] = 1.0 - exp(-1.0*cstar*deltat/tdyn[indsf]*2.0*H2[indsf]/(2.0*H2[indsf] + HI[indsf]))

C  = Ox*Carbon/Oxygen 
Si = Fe*Silicon/Iron

openw,1,filename + 'aux'
FOR i = 0, h.ngas - 1 DO $
   writeu,1,C[i],Ox[i],Si[i],Fe[i];,gsfr[i];,format='($)'
;FOR i = 0, h.ngas - 1 DO $
;   writeu,1,C[i],Ox[i],Si[i],Fe[i],gsfr[i],dummyf[i],dummyf[i],dummyf[i],dummyf[i],dummyi[i]
close,1
stop
END
