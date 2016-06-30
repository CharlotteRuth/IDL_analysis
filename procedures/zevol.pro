PRO zevol,filebase
spawn,'ls ' + filebase + '.param',paramfiles
units = tipsyunits(paramfiles[0])
;spawn,'ls -d ' + filebase + '*.dir',dirnames
readcol,filebase + '.haloid.dat',names,haloid,format='(A,F)'
n = n_elements(dirnames)
n = n_elements(names)
zarray = fltarr(n)
tarray = fltarr(n)
FOR i = 0, n - 1 DO BEGIN
;    cd,dirnames[i]
;    filename = strmid(dirnames[i],0,strlen(dirnames[i]) - 4)
;    print,dirnames[i],', ',filename
    metals = mrdfits(names[i] + '.metals.fits',1)
    ind = where(metals.grp EQ haloid[i])
    print,metals[ind].ox_inneratom,metals[ind].ox_sfr
    rtipsy,names[i],h,g,d,s,/justhead
    tarray[i] = h.time*units.timeunit
    IF metals[ind].ox_sfr GT 0 THEN zarray[i] = metals[ind].ox_sfr ELSE zarray[i] = metals[ind].ox_inneratom
;    cd,'..'
ENDFOR
stop
plot,tarray/1e9,zarray,xtitle = 'Time [Gyr]',ytitle = '12 + log(O/H)',yrange = [6.5,8.5]
IF (keyword_set(outfile)) THEN device,/close ;ELSE stop
END
