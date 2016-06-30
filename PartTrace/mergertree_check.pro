;CC 2/5/14

;This program checks the mergertree by reading in all the virial radii
;and masses and plotting their changes over time.

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

PRO mergertree_check,filebase, finalid = finalid, files = files,halos = halos
IF NOT keyword_set(finalid) THEN finalid = '1'

units = tipsyunits(filebase + '.param')
dist_unit = units.lengthunit
mass_unit = units.massunit
ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)

;find files
IF NOT keyword_set(files) THEN BEGIN
   readcol,filebase + '.grp' + finalid + '.haloid.dat', files, haloids, format='a,l'
   files = reverse(files); + '.amiga.stat'
   haloids = reverse(haloids)
ENDIF
;spawn,'ls */*amiga.stat',files ELSE files = files + '.amiga.stat'
time = fltarr(n_elements(files))
rvir = fltarr(n_elements(files)) 
mvir = fltarr(n_elements(files)) 
FOR i = 0, n_elements(files) - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s,/justhead
    a = h.time
    z = (1 - a)/a
    time[i] = (ageUniverse - wmap3_lookback(z))/1e9 ;in Gyrs
    amigastat = read_stat_struc_amiga(files[i] + '.amiga.stat')
    ind = where(amigastat.group EQ haloids[i])
    rvir[i] = amigastat[ind].rvir
    mvir[i] = amigastat[ind].m_tot
    print,files[i],haloids[i],rvir[i],mvir[i]
 ENDFOR
   set_plot,'x'
   window,0
   plot,time,mvir,psym = -4,xtitle = 'Time',ytitle = 'Virial Mass'
   stop
   plot,time,rvir,psym = -4,xtitle = 'Time',ytitle = 'Virial Radius'
   stop
END
