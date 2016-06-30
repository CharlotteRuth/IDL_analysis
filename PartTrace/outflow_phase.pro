;CC, 12/12/11

;This program determins the characteristics (density, distance from
;center etc.) of the gas at the step
;prior to outflow

;data = fb_gas_character('h516.cosmo25cmb.3072g1MBWK',25000)
FUNCTION outflow_phase,base,totalm = totalm,before = before
loadct,39
spawn,'ls h*param',pfile
units = tipsyunits(pfile)
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
timeunit=SQRT((units.lengthunit*3.086d21)^3/(6.67d-8*units.massunit*1.99d33))/(3600.*24.*365.24)
dens_convert = units.massunit*gm_per_msol*H_per_gm/units.lengthunit^3/cm_per_kpc^3
delta = 1e-3

outflow_z = mrdfits('grp1.rvir.outflow_z.fits')
outflow_iord = mrdfits('grp1.rvir.outflow_iord.fits')
outflow_mass = mrdfits('grp1.rvir.mass_at_outflow.fits')
outflow_disktime = mrdfits('grp1.rvir.outflow_timeindisk.fits')
;stop

;discard all gas particles that are never expelled (z = 99)
;trash = where(outflow_z eq 99.0,complement = keep)
;outflow_z = outflow_z[keep]
;outflow_iord = outflow_iord[keep]
;outflow_mass = outflow_mass[keep]
;outflow_disktime = outflow_disktime[keep]

;Read in inforation about the filenames, halos, time, location and alignment
readcol,'alignment.txt',files,haloid,time,zarray,xc,yc,zc,xa,ya,za,format = '(A,I,F,F,F,F,F,F,F,F)'
center = [[time],[zarray],[xc],[yc],[zc]]
center[*,2] = center[*,2]*units.lengthunit ;1000.0 - units.lengthunit/2.0 ;"center" is in Mpc for a box that goes from [0,0,0] to [units.lengthunit/1000,units.lengthunit/1000,units.lengthunit/1000]
center[*,3] = center[*,3]*units.lengthunit ;*1000.0 - units.lengthunit/2.0
center[*,4] = center[*,4]*units.lengthunit ;*1000.0 - units.lengthunit/2.0
a = [[xa],[ya],[za]]
time0 = 1e9*time[1]/timeunit

;Read in the file names and the halos they correspond to
;readcol,base+'.haloid.dat',files,haloid,format= '(A,F)'
;files = REVERSE(files)          ;going forward in time
;haloid = REVERSE(haloid)

IF keyword_set(before) THEN start = 0 ELSE start = 1
IF keyword_set(before) THEN final = N_ELEMENTS(files) - 2 ELSE final = N_ELEMENTS(files) - 1

;stop

;Iterate through each step, starting from the earliest.  Find out what
;the gas that will be expelled is like at each step
FOR i = 0, N_ELEMENTS(files) - 2 DO BEGIN

;Read in the header for this step and the next step in order to get
;the redshift 
    rtipsy,files[i],h1,g,d,s,/justhead
    scale1 = h1.time
    z1 = 1/scale1 - 1

    rtipsy,files[i + 1],h2,g,d,s,/justhead
    scale2 = h2.time
    z2 = 1/scale2 - 1

;Chose whether we care about the state of the gas before it is
;expelled (earlier step) or after (later step)
    If keyword_set(before) THEN BEGIN
        scale = scale1
        z = z1
    ENDIF ELSE BEGIN
        scale = scale2
        z = z2
    ENDELSE
;Find what gas is lost between this step and the next one
    fb_set2 = where(outflow_z GE z2 - delta AND outflow_z LT z1 - delta)
    fb_set = where(outflow_z GT z2 AND outflow_z NE 99)

;If gas is lost between these two steps continue
    IF fb_set2[0] NE -1 THEN BEGIN
;Select out the data from the outflow information files
        fb_set_iord = outflow_iord[fb_set]
        fb_set_mass = outflow_mass[fb_set]
        fb_set_z    = outflow_z[fb_set]
        fb_set_disktime = outflow_disktime[fb_set]

;Read in the information about the gas at that step (tipsy file, etc)
        IF keyword_set(before) THEN BEGIN
            dir = (strsplit(files[i],'/',/extract))[0]
            title = files[i]
            rtipsy,files[i],h,g,d,s
            undefine,d
;            readarr,files[i] + '.coolontime',h,coolon,part = 'gas',/ascii ;The time at which cooling will turn on again.  I can tell when feedback blast wave last turned off (will turn off).  If I want to know what has experienced feedback between the last two steps, coolon should be since the previous step ( plus a little representing the time period over which cooling was disabled)
;            coolon1 = coolon
;            readarr,files[i + 1] + '.coolontime',h2,coolon2,part = 'gas',/ascii 
            readarr,files[i] + '.iord',h,iord,part = 'gas',/ascii
            iord = long(iord)
            iord1 = iord
            readarr,files[i + 1] + '.iord',h2,iord2,part = 'gas',/ascii 
            iord2 = long(iord2)
        ENDIF ELSE BEGIN
            dir = (strsplit(files[i + 1],'/',/extract))[0]
            title = files[i + 1]
            rtipsy,files[i + 1],h,g,d,s
            undefine,d
;            readarr,files[i + 1] + '.coolontime',h,coolon,part = 'gas',/ascii
;            coolon2 = coolon
;            readarr,files[i] + '.coolontime',h1,coolon1,part = 'gas',/ascii 
            readarr,files[i + 1] + '.iord',h,iord,part = 'gas',/ascii
            iord = long(iord)
            iord2 = iord 
            readarr,files[i] + '.iord',h1,iord1,part = 'gas',/ascii 
            iord1 = long(iord1)
        ENDELSE

;        IF i EQ 0 THEN BEGIN
;            iord_starcoolon = iord1[where(coolon1 NE 0)]
;        ENDIF

        print,dir

;Read in the information about the gas at that step 
        gashistory = mrdfits(dir + '/' + dir + '.allgas.history.fits',1)        
;Match the iords from the gas history file to that of the outflow
        match,fb_set_iord,gashistory.iord,fb_set_ind,gashistory_ind
;Temp is the gashistory information for the outflow
        gfb = gashistory[gashistory_ind]
        fb_set = fb_set[fb_set_ind]
        fb_set_iord = fb_set_iord[fb_set_ind]
        fb_set_mass = fb_set_mass[fb_set_ind]
        fb_set_z    = fb_set_z[fb_set_ind]
        fb_set_disktime = fb_set_disktime[fb_set_ind]

;Match the iords from the outflow to the iords of the tipsy file at
;that output
        match,fb_set_iord,iord,fb_set_ind,iord_ind
;Apply that information to the iord and coolon time array
        iord_fb = iord[iord_ind]
;        coolon_fb = coolon[iord_ind]*timeunit/1e9
        fb_set = fb_set[fb_set_ind]
        fb_set_iord = fb_set_iord[fb_set_ind]
        fb_set_mass = fb_set_mass[fb_set_ind]
        fb_set_z    = fb_set_z[fb_set_ind]
        fb_set_disktime = fb_set_disktime[fb_set_ind]
        gfb = gfb[fb_set_ind]

;Rotations
        IF 1 THEN BEGIN
;Create basis to rotate the particles
            az = reform(a[i,*])
            az1 = a[i,0]
            az2 = a[i,1]
            az3 = a[i,2]
            ax = [az3/sqrt(az1*az1 + az3*az3),0,-1.0*az1/sqrt(az1*az1 + az3*az3)]
            ay = crossp(az,ax) ;[0,-1.0*z/sqrt(y*y + z*z),y/sqrt(y*y + z*z)]
            basis = [[ax],[ay],[az]]

;rotate the outflowing gas from the gashistory file
            gfb.x = (gfb.x - center[i,2])*scale
            gfb.y = (gfb.y - center[i,3])*scale
            gfb.z = (gfb.z - center[i,4])*scale
            gfbpos = [[gfb.x],[gfb.y],[gfb.z]]
            gfbpos = transpose(basis#transpose(gfbpos))
            gfb.x = gfbpos[*,0]
            gfb.y = gfbpos[*,1]
            gfb.z = gfbpos[*,2]

;rotate all the gas particles in the tipsy file
            g.x = (g.x*units.lengthunit - center[i,2])*scale
            g.y = (g.y*units.lengthunit - center[i,3])*scale
            g.z = (g.z*units.lengthunit - center[i,4])*scale
            gpos = [[g.x],[g.y],[g.z]]
            gpos = transpose(basis#transpose(gpos))
            g.x = gpos[*,0]
            g.y = gpos[*,1]
            g.z = gpos[*,2]

;rotate all the stars in the tipsy file
            s.x = (s.x*units.lengthunit - center[i,2])*scale
            s.y = (s.y*units.lengthunit - center[i,3])*scale
            s.z = (s.z*units.lengthunit - center[i,4])*scale
            spos = [[s.x],[s.y],[s.z]]
            spos = transpose(basis#transpose(spos))
            s.x = spos[*,0]
            s.y = spos[*,1]
            s.z = spos[*,2]
        ENDIF ELSE BEGIN
            gfb.x = gfb.x*scale
            gfb.y = gfb.y*scale
            gfb.z = gfb.z*scale
        ENDELSE

;Fix units
;        gfb.vx = gfb.vx*scale
;        gfb.vy = gfb.vy*scale
;        gfb.vz = gfb.vz*scale
        gfb.rho = gfb.rho*dens_convert/scale^3
        g.dens  = g.dens *dens_convert/scale^3

;Add the information about the gas lost at this step to the list of
;gas particles lost
        IF N_ELEMENTS(total) eq 0 THEN BEGIN
            total = gfb
;            totalz = gfbz
;            totalm = gfbm
        ENDIF ELSE BEGIN
            total = [total,gfb]
;            totalz = [totalz,gfbz]
;            totalm = [totalm,gfbm]
        ENDELSE

;Define currently outflowing disk gas (from gas history) that is in the disk
        indcold = where(gfb.temp le 2e4 AND gfb.rho ge 0.1)
;Find the outflow particles (from gas history) that have cooling turned off
;;        indcoolon = where(coolon_fb - time[i] gt 0)
;;        indcoolbefore = where(coolon_fb NE 0)
;       indcoolnext = where(coolon_fb2 NE 0)

;Match the iords of the gas that is outflowing with the tipsy file iords
        match,outflow_iord,iord,outflow_g_ind,g_outflow_ind
;Gas lost selected from the tipsy file
        gas_lost = g[g_outflow_ind]
        gas_lost_iord = iord[g_outflow_ind]
;        gas_lost_coolon = coolon[g_outflow_ind]
;Re-order the disktime array to match the tipsy file
        outflow_disktime_sort = outflow_disktime[outflow_g_ind]
        outflow_z_sort        = outflow_z[outflow_g_ind]
        outflow_iord_sort     = outflow_iord[outflow_g_ind]

;Match the iords of the gas that is _currently_ outflowing with the tipsy file iords
 ;       match,fb_set_iord,iord,fb_g_ind,g_fb_ind
;Gas lost selected from the tipsy file
;        gas_losing = g[g_fb_ind]
;        gas_losing_iord = iord[g_fb_ind]
;        gas_losing_coolon = coolon[g_fb_ind]
;Re-order the disktime array to match the tipsy file
;        fb_set_disktime_sort = fb_set_disktime[outflow_g_ind]
;        fb_set_z_sort        = fb_set_z[outflow_g_ind]

;---------- All gas that is being lost ----------------
;Gas that was last in the disk at this step
        indzdisk =       where(outflow_disktime_sort LT z1 + delta AND outflow_disktime_sort GT z1 - delta)
;Gas that was in the disk for the last time earlier in the sim.
        indzdiskbefore = where(outflow_disktime_sort GT z1 + delta OR outflow_disktime_sort EQ -99 );AND gas_lost_coolon GT 0 AND gas_lost_coolon LE time0)
;Gas that has been or is in the disk
        indzdiskall =    where(outflow_disktime_sort GT z1 - delta)
;Gas that will be in the disk for the last time later on in the simulation
        indzdisklater =  where(outflow_disktime_sort LT z1 - delta)
;Gas that is never in the disk
        indzdisknever =  where(outflow_disktime_sort EQ -99 );AND gas_lost_coolon GT time0 OR gas_lost_coolon EQ 0) 
;        gas_disknever = gas_lost[indzdisknever]
;        iord_disknever = gas_lost_iord[indzdisknever]
;        print,iord_disknever[where(gas_disknever.dens gt 0.1 AND gas_disknever.tempg lt 2e4)

;---------------------- PLOTS -------------
;What particles that are in the disk etc. are currently outflowing?
        match,outflow_iord_sort[indzdisk],      gfb.iord,indzdisk_now,      temp
        match,outflow_iord_sort[indzdisknever], gfb.iord,indzdisknever_now, temp
        match,outflow_iord_sort[indzdiskall],   gfb.iord,indzdiskall_now,   temp
        IF (indzdiskbefore)[0] NE -1 THEN $
          match,outflow_iord_sort[indzdiskbefore],gfb.iord,indzdiskbefore_now,temp
        IF (indzdisklater)[0]  NE -1 THEN $
          match,outflow_iord_sort[indzdisklater], gfb.iord,indzdisklater_now, temp

        IF 1 THEN BEGIN
;How many of the gas particles that have never been in the disk and are
;about to be expelled have been supernova heated?
;            match,outflow_iord_sort[indzdisknever[indzdisknever_now]],gfb[indcoolbefore].iord,indzdisknever_now_coolon,temp
;            print,'Fraction of outflowing gas never in the disk that were heated by SN: ' + $
;              STRTRIM(N_ELEMENTS(indzdisknever_now_coolon),2) + ' / ' + strtrim(N_ELEMENTS(indzdisknever_now),2) + ' = ' + strtrim(double(n_elements(indzdisknever_now_coolon))/n_elements(indzdisknever_now),2)
;And the fraction of gas that has been in the disk
;            IF (indzdiskall)[0] NE -1  THEN BEGIN
;                match,outflow_iord_sort[indzdiskall[indzdiskall_now]],gfb[indcoolbefore].iord,indzdiskall_now_coolon,temp
;                print,'Fraction of outflowing gas once in the disk that were heated by SN: ' + $
;                  STRTRIM(N_ELEMENTS(indzdiskall_now_coolon),2) + ' / ' + strtrim(N_ELEMENTS(indzdiskall_now),2) + ' = ' + strtrim(double(n_elements(indzdiskall_now_coolon))/n_elements(indzdiskall_now),2)
;            ENDIF
;            IF (indzdiskall)[0] NE -1  THEN BEGIN
;                match,outflow_iord_sort[indzdiskall[indzdiskall_now]],gfb[indcoolbefore].iord,indzdiskall_now_coolon,temp
;                print,'Fraction of outflowing gas that was once in the disk: ' + $
;                  strtrim(n_elements(indzdiskall_now),2) + ' / ' + strtrim(n_elements(fb_set),2) + ' = ' + strtrim(double(n_elements(indzdiskall_now))/n_elements(fb_set),2)
;            ENDIF
;            match,outflow_iord_sort[indzdiskall[indzdiskall_now]],gfb[indcoolbefore].iord,indzdiskall_now_coolon,temp
;            print,'Fraction of outflowing gas that was heated by supernovae: ' + $
;              strtrim(n_elements(where(coolon_fb NE 0)),2) + ' / ' + strtrim(n_elements(fb_set),2) + ' = ' + strtrim(double(n_elements(where(coolon_fb NE 0)))/n_elements(fb_set),2)
;And the gas that is currently in the disk
        ENDIF

        IF 0 THEN BEGIN
;Do I need to worry about gas that was once in the disk and has been
;expelled but has never been heated by a supernova?
;Distribution of when the outflowing gas was last in the disk
            window,4
            histogramp,outflow_disktime_sort[indzdiskall[indzdiskall_now]],min = 1e-6,max = 5,nbins = 25,xtilte =' Redshift',ytitle = 'Amount of gas lost from disk'
            histogramp,outflow_disktime_sort[indzdiskall[indzdiskall_now[indzdiskall_now_coolon]]],min = 1e-6,max = 5,nbins = 25,/overplot,linestyle = 2
;What I find is that for the gas that had been in the disk and now is
;gone: If it wasn't heated by SN, it was more likely to have been in
;the disk a long time ago

;Probably what I really want to know is: how much gas has been lost from the
;disk recently? 
        ENDIF

        IF 1 THEN BEGIN
            window,0,xsize = 600,ysize = 600
            plot,gfb.x,gfb.z,psym = 3,xtitle = 'X [kpc]',ytitle = 'Z [kpc]',/nodata,title = title
            oplot,s.x,s.z,psym = 3
            oplot,gfb.x,           gfb.z,           color = 210,psym = 3
            IF (indcold)[0] NE -1 THEN oplot,gfb[indcold].x,  gfb[indcold].z,  color = 190,psym = 3
;            IF (indcoolon)[0] NE -1 THEN $
;              oplot,gfb[indcoolon].x,gfb[indcoolon].z,color = 240,psym = 3
            oplot,gas_lost[indzdisknever[indzdisknever_now]].x,   gas_lost[indzdisknever[indzdisknever_now]].z, color = 150,psym = 3
            if (indzdiskbefore)[0] NE -1 THEN $
              IF (indzdiskbefore_now)[0] NE -1 THEN $
              oplot,gas_lost[indzdiskbefore[indzdiskbefore_now]].x,gas_lost[indzdiskbefore[indzdiskbefore_now]].z,color = 60,psym = 3
            IF (indzdisklater)[0]  NE -1 THEN $
              IF (indzdisklater_now)[0]  NE -1 THEN $
              oplot,gas_lost[indzdisklater[indzdisklater_now]].x, gas_lost[indzdisklater[indzdisklater_now]].z, color = 100, psym = 3
            IF (indzdisk_now)[0]  NE -1 THEN $
              oplot,gas_lost[indzdisk[indzdisk_now]].x,        gas_lost[indzdisk[indzdisk_now]].z,      color = 60, psym = 3
;            stop
;            oplot,gfb[indcoolbefore].x,gfb[indcoolbefore].z,color = 240,psym = 3
        ENDIF

        IF 1 THEN BEGIN
            window,1
            plot,gfb.rho,            gfb.temp,title = title, xtitle = 'Density [amu/cc]',ytitle = 'Temp [K]',/xlog,/ylog,yrange = [1e2,1e8],xrange = [1e-6,1e4],/nodata
            oplot,gfb.rho,           gfb.temp,           psym = 3,color = 210
            IF (indcold)[0] NE -1 THEN oplot,gfb[indcold].rho,  gfb[indcold].temp,  psym = 3,color = 190
;            IF (indcoolon)[0] NE -1 THEN $
;              oplot,gfb[indcoolon].rho,gfb[indcoolon].temp,psym = 3,color = 240
;Gas that will never be in the disk
            oplot,gas_lost[indzdisknever[indzdisknever_now]].dens,   gas_lost[indzdisknever[indzdisknever_now]].tempg, color = 150,psym = 3 
;Gas that was in the disk for the last time earlier on
            IF (indzdiskbefore)[0] NE -1 THEN $
              IF (indzdiskbefore_now)[0] NE -1 THEN $
              oplot,gas_lost[indzdiskbefore[indzdiskbefore_now]].dens,gas_lost[indzdiskbefore[indzdiskbefore_now]].tempg,color = 60,psym = 3 
;Gas that will be in the disk for the last time later
            IF (indzdisklater)[0]  NE -1 THEN $
              If (indzdisklater_now)[0]  NE -1 THEN $
              oplot,gas_lost[indzdisklater[indzdisklater_now]].dens, gas_lost[indzdisklater[indzdisklater_now]].tempg, color = 100,psym = 3 
;Gas that is currently in the disk for the last time before being expelled
            IF (indzdisk_now)[0]  NE -1 THEN $
              oplot,gas_lost[indzdisk[indzdisk_now]].dens,        gas_lost[indzdisk[indzdisk_now]].tempg,      color = 60,psym = 3 
;            stop
;            oplot,gfb[indcoolbefore].rho,gfb[indcoolbefore].temp,psym = 3,color = 240
        ENDIF
    
        IF 0 THEN BEGIN
            window,2
            indz = where(outflow_z lt z2 + delta AND outflow_z GT z2 - delta)
            histogramp,outflow_disktime[indz],min = -0.001, max = 6, nbins = 100, xrange = [0,6],title = 'Redshift When Current Outflow Last in Disk'
        ENDIF

        IF 0 THEN BEGIN
            window,3
            histogramp,outflow_z_sort[indzdisknever],min = -0.001,max = 6, nbins = 100, xrange = [0,6],title = 'Redshift When Gas Will Exit Disk'
            histogramp,outflow_z_sort[indzdisknever],   /overplot,color = 150,min = -0.001,max = 6, nbins = 100
            IF (indzdiskbefore)[0] NE -1 THEN $
              histogramp,outflow_z_sort[indzdiskbefore],/overplot,color = 100,min = -0.001,max = 6, nbins = 100
            IF (indzdisklater)[0]  NE -1 THEN $
              histogramp,outflow_z_sort[indzdisklater], /overplot,color =  30,min = -0.001,max = 6, nbins = 100
            histogramp,outflow_z_sort[indzdisk],        /overplot,color =  60,min = -0.001,max = 6, nbins = 100
        ENDIF
        
        stop

    ENDIF
ENDFOR

return,total
END

PRO outflow_phase_master,outplot = outplot
files = ['h516.cosmo25cmb.3072g1MBWK','h516.cosmo25cmb.3072g14HBWK']
dir = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK','/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK']
kpcunit = [25000.,25000.]
colors = [50,240]

files = ['h986.cosmo50cmb.3072g1bwK','h986.cosmo50cmb.3072gs1MbwK','h986.cosmo50cmb.3072g14HBWK']
dir = ['/astro/store/nbody3/fabio/h986/3072g1bwK',$
       '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK',$
       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK']
kpcunit = [50000.,50000.,50000.]
msolunit = [1.84793e16,1.84793e16,1.84793e16]
colors = [50,120,240]
thicks = [6,6,6]

;files = ['h603.cosmo50cmb.3072g1bwK','h603.cosmo50cmb.3072gs1MbwK','h603.cosmo50cmb.3072g14HBWK']
;dir = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.3072g1bwK',$
;       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK',$
;       '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK']
;kpcunit = [50000.,50000.,50000.]
;Msolunit = [1.84793e16,1.84793e16,1.84793e16]
;colors = [50,120,240]

formatplot,outplot=outplot,/thick
IF KEYWORD_SET(outplot) THEN device,filename = outplot,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch,/color

FOR i = 2, N_ELEMENTS(files)-1 DO BEGIN
    cd,dir[i]
    fb_gas = outflow_phase(files[i],kpcunit[i],msolunit[i],before = 0)
    IF i EQ 0 THEN plot,fb_gas.rho,fb_gas.temp,xtitle = 'Density [amu/cc]',ytitle = 'Temp [K]',/xlog,/ylog,/nodata,yrange = [1e2,1e8],xrange = [1e-6,1e4]
;    oplot,fb_gas.rho,fb_gas.temp,psym = 3,color = colors[i]
;    histogramp,SQRT(fb_gas.x*fb_gas.x + fb_gas.y*fb_gas.y),/overplot,/normalize,color = colors[i],thick = thicks[i]
    stop
ENDFOR

IF KEYWORD_SET(outplot) THEN device,/close

END
