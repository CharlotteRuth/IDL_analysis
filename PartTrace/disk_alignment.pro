;CC 12/13/11

;This program finds the alignment of the disk of the galaxy for each
;timestep and outputs to a file: 
;(file, haloid, t, z, x_c, y_c, z_c,a_x, a_y, a_z)
; where <a_x,a_y,a_z> is the vector defining the plane of the galaxy
;.r tipsySatsHI
;.r alignHI
;readcol, base + '.haloid.dat', files, haloid, format= '(A,F)'
;files = reverse(files) ;going forward in time
;haloid = reverse(haloid)
;data = disk_alignment(filebase, finalid = finalid, files = files, halos = haloid)

PRO disk_alignment,filebase, finalid = finalid, files = files,halos = halos
IF NOT keyword_set(finalid) THEN finalid = '1'

units = tipsyunits(filebase + '.param')
dist_unit = units.lengthunit
mass_unit = units.massunit
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)

;find files
IF NOT keyword_set(files) THEN BEGIN
   readcol,filebase + '.grp' + finalid + '.haloid.dat', files, halos, format='a,l'
   files = reverse(files); + '.amiga.stat'
   halos = reverse(halos)
ENDIF
;spawn,'ls */*amiga.stat',files ELSE files = files + '.amiga.stat'
center = fltarr(n_elements(files),11)
close,11
openw,11,'grp' + finalid + '.alignment.txt'
FOR i = 0, n_elements(files) - 1 DO BEGIN
    print,files[i]
;    endpos = strpos(files[i],'/')
;    base = strmid(files[i],0,endpos)
;    dotpos = strsplit(files[i],'.')
;    base = strmid(files[i],0,dotpos[n_elements(dotpos) - 2] - 1)
    base = files[i]
    rtipsy,files[i],h,g,d,s,/justhead
    a = h.time
    z = (1 - a)/a
    time = (ageUniverse - wmap3_lookback(z))/1e9 ;in Gyrs
    satsdata = read_stat_struc_AMIGA(files[i] + '.amiga.stat')
    IF keyword_set(halos) THEN j = where(satsdata.group EQ halos[i]) ELSE j = 0
    
    tipsysatshi, files[i], fix(halos[i]), dist_unit, mass_unit,cutout_rad = 1,notipsy = 1,spinaxes = spinaxes,cx = cx, cy = cy, cz = cz,valign=valign
;    cx = cx*dist_unit/1000.0
;    cy = cy*dist_unit/1000.0
;    cz = cz*dist_unit/1000.0
    center[i,*] = [time,z,cx,cy,cz,valign,spinaxes]
;    print,base,fix(halos[i]),time,z,cx,cy,cz,valign,spinaxes,format = '(A,I,F,F,F,F,F,F,F,F)'
    printf,11,base,fix(halos[i]),time,z,cx,cy,cz,valign,spinaxes,format = '(A,I,F,F,F,F,F,F,F,F,F,F,F)'
ENDFOR
close,11
;RETURN,center
END
