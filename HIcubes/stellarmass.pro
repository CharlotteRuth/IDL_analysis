;cd,'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00072.dir'
;file = 'h603.cosmo50cmb.3072g14HBWK.00072'
;pfile = '../../h603.cosmo50cmb.3072g14HBWK.param'
;nhalos = 3265


pro stellarmass,file,pfile,nhalos = nhalos,mint = mint
rtipsy,file,h,g,d,s
read_tipsy_arr,file + '.amiga.grp',h,stargrp,part = 'star'
units = tipsyunits(pfile)
IF NOT KEYWORD_SET(mint) THEN mint = 0.6*1e9 ;Gyr

FMT = 'I,X,X,X,X,F,X,F,F,F,F,X,F,X,X,X,X,X,X,A,A'
readcol,file + '.amiga.stat',F=FMT,halon,vir,gas,star,dark,velmax,veldisp,cont,sat,/silent
IF ((where(cont eq 'clean'))[0] eq -1) THEN BEGIN
    FMT = 'I,X,X,X,X,F,X,F,F,F,F,X,F,X,X,X,X,X,X,X,A,A'
    readcol,file + '.amiga.stat',F=FMT,halon,vir,gas,star,dark,velmax,veldisp,cont,sat,/silent
ENDIF
IF ((where(sat eq 'no'))[0] eq -1) THEN BEGIN
    FMT = 'I,X,X,X,X,F,X,F,F,F,F,X,F,X,X,X,X,X,X,A'
    readcol,file + '.amiga.stat',F=FMT,halon,vir,gas,star,dark,velmax,veldisp,cont,/silent
ENDIF
stop
IF not keyword_set(nhalos) THEN nhalos = n_elements(halon)
ystar = fltarr(nhalos)

FOR i = 0, nhalos - 1 DO BEGIN
    ind = where(stargrp eq halon[i] AND s.tform*units.timeunit gt mint)
    IF ind[0] ne -1 THEN BEGIN
        ystar[i] = TOTAL(s[ind].mass)*units.massunit
    ENDIF Else ystar[i] = 0
ENDFOR
stop
writecol,'h603.cosmo50cmb.3072g14HBWK.00072.halodat',halon,vir,dark,gas,star,ystar
end


;FOR i = 0, nhalos - 1 DO IF (where(stargrp eq halon[i] AND s.tform*units.timeunit gt mint))[0] ne -1 THEN ystar[i] = TOTAL(s[where(stargrp eq halon[i] AND s.tform*units.timeunit gt mint)].mass)*units.massunit ELSE ystar[i] = 0
