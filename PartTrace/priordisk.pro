;This procedure finds when particles were last in the disk before
;being expelled by SN

PRO priordisk
loadct,39

ageUniverse = 13.7346*1e9 ;wmap3_lookback(10000)
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

delta = 1e-5
mindens = 0.1 ;amu/cc minmum density of gas defined to be in the disk
maxtemp = 2e4; maximum temperatures of gas defined to be in the disk

spawn,'ls h*param',pfile
units = tipsyunits(pfile)
timeunit=SQRT((units.lengthunit*3.086d21)^3/(6.67d-8*units.massunit*1.99d33))/(3600.*24.*365.24)
dens_convert = units.massunit*gm_per_msol*H_per_gm/units.lengthunit^3/cm_per_kpc^3
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

;Select data that is currently outflowing
    outflow_ind = where(outflow_z GE znext - delta AND outflow_z LT z - delta)
;If there is gas that is currently outflow, look at it
    IF outflow_ind[0] ne -1 THEN BEGIN
        IF keyword_set(ejected_iord) THEN match,ejected_iord,outflow_iord[outflow_ind],subejet,suboutflow
        IF NOT keyword_set(ejected_iord) THEN $
          ejected_iord = outflow_iord[outflow_ind] ELSE $
          ejected_iord = [ejected_iord,outflow_iord[outflow_ind]]
; Add onto the list of particles that have been ejected all those that
; were ejected between this current step and the one after it (in time)

        g = mrdfits(files[i] + '.allgas.history.fits',1)
;Calculate the center of the halo of interest
        stat = read_stat_struc_amiga(files[i] + '.amiga.stat')
        ind = where(stat.group eq haloid[i])
        stat = stat[ind]
        stat.xc = stat.xc*1000.0 - units.lengthunit/2
        stat.yc = stat.yc*1000.0 - units.lengthunit/2
        stat.zc = stat.zc*1000.0 - units.lengthunit/2
        g.rho = g.rho*dens_convert/h.time^3
        g.x = g.x - stat.xc
        g.y = g.y - stat.yc
        g.z = g.z - stat.zc
        gall = g
        
        match,ejected_iord,g.iord,ejected_iord_ind,iord_ind
        IF n_elements(ejected_iord) NE n_elements(ejected_iord_ind) THEN BEGIN
            print,'problem'
            stop
        ENDIF
        ejected_iord = ejected_iord[ejected_iord_ind]
        g = g[iord_ind]
        disk_ind = where(g.rho GE mindens AND g.temp LE maxtemp AND sqrt(g.x^2 + g.y^2 + g.z^2) LE stat.rvir, complement = nodisk_ind)

;        rtipsy,files[i],h0,g0,d0,s0
;        readarr,files[i] + '.iord',h0,iord,part = 'gas',/ascii 
;        g0.dens = g0.dens*dens_convert/h.time^3
;        g0.x = g0.x*units.lengthunit - stat.xc
;        g0.y = g0.y*units.lengthunit - stat.yc
;        g0.z = g0.z*units.lengthunit - stat.zc
        
;        match,ejected_iord,iord,ejected_iord_ind0,iord_ind0
;        ejected_iord0 = ejected_iord[ejected_iord_ind0]
;        g0 = g0[iord_ind0]       
;        disk_ind0 = where(g0.dens GE mindens AND g0.tempg LE maxtemp AND sqrt(g0.x^2 + g.y^2 + g.z^2) LE stat.rvir, complement = nodisk_ind0)        

        IF 0 THEN BEGIN
            window,0,xsize = 600,ysize = 600
            plot,g.x,g.z,psym = 3,xtitle = 'X [kpc]',ytitle = 'Y [kpc]',title = files[i]
            oplot,g[nodisk_ind].x,g[nodisk_ind].z,psym = 3,color = 240
            IF (disk_ind)[0] ne -1 THEN oplot,g[disk_ind].x,g[disk_ind].z,psym = 3,color = 60
;            oplot,g0[nodisk_ind0].x,g0[nodisk_ind0].z,psym = 3,color = 200
;            IF (disk_ind0)[0] ne -1 THEN oplot,g0[disk_ind0].x,g0[disk_ind0].z,psym = 3,color = 100
            
            window,1
            plot,g.rho,g.temp,psym = 3,/xlog,/ylog,title = files[i],xrange = [1e-6,1e4],yrange = [100,1e6] 
;            oplot,gall.rho,gall.temp,psym = 3
            oplot,g[nodisk_ind].rho,g[nodisk_ind].temp,psym = 3,color = 240
            IF (disk_ind)[0] ne -1 THEN oplot,g[disk_ind].rho,g[disk_ind].temp,psym = 3,color = 60
;            oplot,g0[nodisk_ind0].dens,g0[nodisk_ind0].tempg,psym = 3,color = 200
;            IF (disk_ind0)[0] ne -1 THEN oplot,g0[disk_ind0].dens,g0[disk_ind0].tempg,psym = 3,color = 100
        ENDIF
            
        IF (disk_ind)[0] ne -1 THEN BEGIN
            match,outflow_iord_base,g[disk_ind].iord,  outflow_iord_ind,temp
            outflow_timeindisk_base[outflow_iord_ind] = z
        ENDIF

        print,N_ELEMENTS(ejected_iord),N_ELEMENTS(nodisk_ind), N_ELEMENTS(disk_ind),N_ELEMENTS(where(outflow_timeindisk_base ne 0))
        match,ejected_iord,g[nodisk_ind].iord,ejected_iord_ind,temp
        ejected_iord = ejected_iord[ejected_iord_ind]
        stop
    ENDIF
    print,files[i]
ENDFOR
match,outflow_iord_base,ejected_iord,outflow_iord_ind,temp
outflow_timeindisk_base[outflow_iord_ind] = -99
stop
mwrfits,outflow_timeindisk_base,'grp1.rvir.outflow_timeindisk.fits',/create

END
