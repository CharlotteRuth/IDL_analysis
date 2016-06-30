;This procedure finds when particles were last in the disk before
;being expelled by SN

PRO priordisk,base,kpcunit,msolunit
loadct,39

ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
;timeunit=SQRT((kpcunit*3.086d21)^3/(6.67d-8*msolunit*1.99d33))/(3600.*24.*365.24)
dens_convert = msolunit*gm_per_msol*H_per_gm/(kpcunit)^3/cm_per_kpc^3;/scale^3

delta = 1e-5
mindens = 0.1 ;amu/cc minmum density of gas defined to be in the disk
maxtemp = 2e4; maximum temperatures of gas defined to be in the disk

outflow_z_base = mrdfits('grp1.rvir.outflow_z.fits')
outflow_iord_base = mrdfits('grp1.rvir.outflow_iord.fits')
outflow_timeindisk_base = outflow_z_base*0 ;Last time at which ejected gas is in the disk
trash = where(outflow_z_base eq 99.0,complement = keep) ;Get rid of all particles that are never ejected
outflow_z = outflow_z_base[keep] ;The redshift at which the gas is ejected
outflow_iord = outflow_iord_base[keep] ;The iord values of the ejected gas

;readcol,base+'.haloid.dat',files,haloid,format= '(A,F)'
readcol,'alignment.txt',files,haloid,time,z,xc,yc,zc,xa,ya,za,format = '(A,I,F,F,F,F,F,F,F,F)'
files = REVERSE(files)          ;going forward in time
haloid = REVERSE(haloid)

FOR i = 1, N_ELEMENTS(files) - 1 DO BEGIN
;Going backwards in time
    rtipsy,files[i - 1],hnext,g,d,s,/justhead
    znext = 1/hnext.time - 1

    rtipsy,files[i],h,g,d,s,/justhead
    z = 1/h.time - 1

    outflow_ind = where(outflow_z GE znext - delta AND outflow_z LT z - delta)
    IF outflow_ind[0] ne -1 THEN BEGIN
        IF keyword_set(ejected_iord) THEN match,ejected_iord,outflow_iord[outflow_ind],subejet,suboutflow
 ;       stop
        IF NOT keyword_set(ejected_iord) THEN $
          ejected_iord = outflow_iord[outflow_ind] ELSE $
          ejected_iord = [ejected_iord,outflow_iord[outflow_ind]]
; Add onto the list of particles that have been ejected all those that
; were ejected between this current step and the one after it (in time)

        rtipsy,files[i],h,g,d,s
        readarr,files[i] + '.iord',h,iord,part = 'gas',/ascii
        readarr,files[i] + '.amiga.grp',h,grp,part = 'gas',/ascii
        
  ;      match2,ejected_iord,iord,ejected_iord_ind,iord_ind ;iord_ind = the indicies of g, iord, and grp that correspond to the ejected particles
  ;      iord_ind = iord_ind[where(iord_ind ne -1)]
        match,ejected_iord,iord,ejected_iord_ind,iord_ind

        IF n_elements(ejected_iord) NE n_elements(ejected_iord_ind) THEN BEGIN
            print,'problem'
            stop
        ENDIF
        ejected_iord = ejected_iord[ejected_iord_ind]
        g = g[iord_ind]
        iord = iord[iord_ind]
        grp = grp[iord_ind]

;       disk_ind = where(grp eq haloid[i] AND g.dens*dens_convert/h.time^3 ge mindens AND g.tempg le maxtemp, complement = noejec_ind) ; select all gas that is currently in the disk (dense enough and cold enough) of the primary galaxy
        disk_ind = where(g.dens*dens_convert/h.time^3 ge mindens AND g.tempg le maxtemp, complement = nodisk_ind)
        
        IF 0 THEN BEGIN
            window,0,xsize = 600,ysize = 600
            plot,g.x*kpcunit*h.time,g.z*kpcunit*h.time,psym = 3,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',title = files[i]
            oplot,g[nodisk_ind].x*kpcunit*h.time,g[nodisk_ind].z*kpcunit*h.time,psym = 3,color = 240
            IF (disk_ind)[0] ne -1 THEN oplot,g[disk_ind].x*kpcunit*h.time,g[disk_ind].z*kpcunit*h.time,psym = 3,color = 60
            
            window,1
            plot,g.dens*dens_convert/h.time^3,g.tempg,psym = 3,/xlog,/ylog,title = files[i] 
            oplot,g[nodisk_ind].dens*dens_convert/h.time^3,g[nodisk_ind].tempg,psym = 3,color = 240
            IF (disk_ind)[0] ne -1 THEN oplot,g[disk_ind].dens*dens_convert/h.time^3,g[disk_ind].tempg,psym = 3,color = 60
        ENDIF
            
        IF (disk_ind)[0] ne -1 THEN BEGIN
            match,outflow_iord_base,iord[disk_ind],  outflow_iord_ind,temp
            outflow_timeindisk_base[outflow_iord_ind] = z
        ENDIF

;        help,ejected_iord
        match,ejected_iord,     iord[nodisk_ind],ejected_iord_ind,temp
        ejected_iord = ejected_iord[ejected_iord_ind]
        print,files[i]
;        help,disk_ind
;        help,ejected_iord
;        stop
    ENDIF
;    stop
ENDFOR

match,outflow_iord_base,ejected_iord,outflow_iord_ind,temp
outflow_timeindisk_base[outflow_iord_ind] = -99

mwrfits,outflow_timeindisk_base,'grp1.rvir.outflow_timeindisk.fits',/create

;stop
END
