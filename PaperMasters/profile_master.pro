;profile_master,outplot='~/Plots/profile/'
PRO profile_master,outplot = outplot,color = color,verbose = verbose,formatthick = formatthick,debug = debug,galtotal = galtotal
formatplot,outplot = outplot,thick = formatthick
IF keyword_set(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF keyword_set(outplot) THEN bgcolor = 255 ELSE bgcolor = 0

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/nobackupp2/crchrist/MolecH/' $
ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
IF hostname EQ 'ozma' THEN datafile_prefix = '/home/christensen/Code/Datafiles/' ELSE datafile_prefix = '/astro/users/christec/code/Datafiles/'
zsolar = 0.0130215

pre_cuthalos = 0              ;pre
plot_profile = 0
do_metal = 0
write_metal = 0
plot_schmidtlaw_global_obs = 0
halo_info = 0
sfh = 1

mu_c50 = 1.84793e16
mu_c25 = 2.310e15
lu_c50 = 50000.
lu_c25 = 25000.

x = 1
CASE x OF
   1: BEGIN
        dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
        file799 = 'h799.cosmo25cmb.3072g14HBWK'
        key799 = 'h799'
        dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
        file516 = 'h516.cosmo25cmb.3072g14HBWK'
        key516 = 'h516'
        dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
        file986 = 'h986.cosmo50cmb.3072g14HBWK'
        key986 = 'h986'
        dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
        file603 = 'h603.cosmo50cmb.3072g14HBWK'
        key603 = 'h603'
        dir277 = prefix + 'h277.cosmo50cmb.3072g/h277.cosmo50cmb.3072g14HMbwK/'
        file277 = 'h277.cosmo50cmb.3072g14HMbwK'
        key277 = 'h277'
        dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
        file258 = 'h258.cosmo50cmb.3072g14HMbwK'
        key258 = 'h258'
        dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
        file285 = 'h285.cosmo50cmb.3072g14HMbwK'
        key285 = 'h285'
        dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
        file239 = 'h239.cosmo50cmb.3072g14HMbwK'
        key239 = 'h239'
        dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277,dir277,dir258,dir258,dir285,dir285,dir285,dir239]
        files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277,file277,file258,file258,file285,file285,file285,file239]
        steps = [['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00032','0'    ,'0'    ,'00080','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00032','0'    ,'0'    ,'00080','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['0'    ,'0'    ,'0'    ,'0'    ,'00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['0'    ,'0'    ,'0'    ,'0'    ,'00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['0'    ,'0'    ,'0'    ,'0'    ,'00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512'],$
                 ['00024','00036','00048','00060','00084','00132','00168','00228','00276','00324','00408','00512']]
;        step = '00512'
;        outfiles = dirs + files + '.' + step + '/' + files + '.' + step
        nz = (size(steps))[1]
        ngal = (size(steps))[2]
        haloid = ['1'   ,$ ;h799
                  '1'   ,'2'   ,$ ;h516
                  '1'   ,'2'   ,'3'  ,'8'   ,$ ;h986
                  '1'   ,'2'   ,'3'  ,$ ;h603
                  '1'   ,'2'   ,$ ;h277
                  '1'   ,'4'   ,$ ;h258
                  '1'   ,'4'   ,'9'  ,$ ;h285
                  '1'] ;h239
        dwarfs = [0, 1, 2, 4, 5, 6, 10, 11, 13, 15, 16]
        distunits = [fltarr(3) + lu_c25, fltarr(15) + lu_c50]
        massunits = [fltarr(3) + mu_c25, fltarr(15) + mu_c50]
        key  = [key799,key516,key516,key986,key986,key986,key986,key603,key603,key603,key277,key277,key258,key258,key285,key285,key285,key239] + ', ' + haloid
;        dirs    = [dir799,dir516,dir516,dir986,dir986,dir986,dir986,dir603,dir603,dir603,dir277]
;        files    = [file799,file516,file516,file986,file986,file986,file986,file603,file603,file603,file277]
;        haloid = ['1'   ,'1'   ,'2'   ,'1'   ,'2'   ,'3'   ,'9'   ,'1',   '2'    ,'3'  ,'2']

;        colors = fltarr(n_elements(files)) + fgcolor
;        loadct = fltarr(n_elements(files))
;        linestyles = fltarr(n_elements(files))
        psym = fltarr(n_elements(files)) + 16
        psym[dwarfs] = 17
        obscolor = 100
        obssym = 2
    END
ENDCASE

filebases = [(files + '.halo.' + haloid)]

n = n_elements(files)
IF keyword_set(color) THEN BEGIN
    loadct,39
    IF NOT keyword_set(ctables) THEN ctables = [39,39,39]
    IF NOT keyword_set(obscolor) THEN obscolor = fgcolor
    IF NOT keyword_set(colors) THEN  colors  = (findgen(n) + 1)*254/n ELSE colors = colors
    IF NOT keyword_set(psym) THEN psym = fltarr(n) + 4
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2
    IF NOT keyword_set(linestyles) THEN linestyles = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 1.5
ENDIF ELSE BEGIN
    loadct,0    
    IF NOT keyword_set(ctables) THEN ctables = [0,0,0]
    IF NOT keyword_set(obscolor) THEN obscolor = 100
    If NOT keyword_set(colors) THEN  colors = (fltarr(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT keyword_set(psym) THEN  psym = (findgen(n)+2)*2
    IF NOT keyword_set(thicks) THEN thicks = fltarr(n) + 2 ;thicks = (findgen(n) + 1)*6/n - 1
    IF NOT keyword_set(linestyles) THEN linestyles = REVERSE(findgen(n)*2) 
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 1.5
ENDELSE

IF keyword_set(outplot) THEN BEGIN
    xsize = 18                  ;10*n
    ysize = 12                  ;12
ENDIF ELSE BEGIN
    xsize = 400                 ;*n
    ysize = 350                 ;475
ENDELSE

; ------------------------------------------ Cut out halo ----------------------------
IF pre_cuthalos THEN BEGIN
    outarray =['amiga.grp','FeMassFrac','OxMassFrac','HI','H2','coolontime','iord'];,'HeI','HeII']
    FOR i = 0, n_elements(dirs) - 1 DO BEGIN
        valid = where(steps[*,i] NE '0')
        FOR j = 0, nz - 1 DO BEGIN
            IF (where(j EQ valid))[0] NE -1 THEN BEGIN
                readcol,dirs[i] + files[i] + '.grp' + haloid[i] + '.haloid.dat',stepnames,stephaloids,format='(A,A)'
                k = where(stepnames EQ files[i] + '.' + steps[j,i] + '/' + files[i] + '.' + steps[j,i])
                outfile = dirs[i] + files[i] + '.' + steps[j,i] + '/' + files[i] + '.' + steps[j,i]
                IF NOT (file_test(outfile + '.halo.' + stephaloids[k] + '.std') AND $
                        file_test(outfile + '.halo.' + stephaloids[k] + '.amiga.grp') AND $
                        file_test(outfile + '.halo.' + stephaloids[k] + '.FeMassFrac') AND $
                        file_test(outfile + '.halo.' + stephaloids[k] + '.OxMassFrac') AND $
                        file_test(outfile + '.halo.' + stephaloids[k] + '.HI') AND $
                        file_test(outfile + '.halo.' + stephaloids[k] + '.H2') AND $
                        file_test(outfile + '.halo.' + stephaloids[k] + '.coolontime') AND $
                        file_test(outfile + '.halo.' + stephaloids[k] + '.iord')) THEN BEGIN
; AND $
;                        file_test(outfile + '.halo.' + stephaloids[k] + '.HeI') AND $
;                        file_test(outfile + '.halo.' + stephaloids[k] + '.HeII')) THEN BEGIN
                    print,outfile + '.halo.' + stephaloids[k]
                    tipsysatshi, outfile, stephaloids[k], distunits[i], massunits[i], /cutout_rad, outarray = outarray,/gethi
                    close,/all
                ENDIF
            ENDIF
        ENDFOR
    ENDFOR
ENDIF

;-------------------------------------------- plot_profile ------------------------------
; FOR i = 0, ngal - 1 DO BEGIN
IF plot_profile THEN BEGIN
    FOR i = 0, n_elements(dirs) - 1 DO BEGIN
;    FOR i = 17, 17 DO BEGIN
        valid = where(steps[*,i] NE '0')
        validids = strarr(n_elements(valid))
        FOR j = 0, n_elements(valid) - 1 DO BEGIN
;        FOR j = 10, n_elements(valid) - 1 DO BEGIN
            readcol,dirs[i] + files[i] + '.grp' + haloid[i] + '.haloid.dat',stepnames,stephaloids,format='(A,A)'
            k = where(stepnames EQ files[i] + '.' + steps[valid[j],i] + '/' + files[i] + '.' + steps[valid[j],i])            
            validids[j] = stephaloids[k]
        ENDFOR
        ; rather than haloids[i], use validids[j]
        IF keyword_set(outplot) THEN plot_profile,$
          dirs[i] + files[i] + '.' + steps[valid,i]+ '/' + files[i] + '.' + steps[valid,i] + '.halo.' + validids,$
          fltarr(nz) + massunits,$
          fltarr(nz) + distunits,$
          filebase = files[i] + '.' + steps[valid,i] + '.halo.' + haloid[i],$
          outplot = outplot + '/' + files[i] + '.' + steps[valid,i] + '.halo.' + haloid[i], $
          debug = debug,/starlog,galtotal = galtotal $
        ELSE plot_profile,$
          dirs[i] + files[i] + '.' + steps[valid,i] + '/' + files[i] + '.' + steps[valid,i] + '.halo.' + validids,$
          fltarr(nz) + massunits,$
          fltarr(nz) + distunits,$
          filebase = files[i] + '.' + steps[valid,i] + '.halo.' + haloid[i],$
          debug = debug,/starlog
    ENDFOR
ENDIF

;------------------------------------------ halo_info ----------------------------------
IF halo_info THEN BEGIN
    FOR i = 1, n_elements(dirs) - 1 DO BEGIN
;    FOR i = 17, 17 DO BEGIN
        valid = where(steps[*,i] NE '0')
        validids = strarr(n_elements(valid))
        FOR j = 0, n_elements(valid) - 1 DO BEGIN
;        FOR j = 10, n_elements(valid) - 1 DO BEGIN
            readcol,dirs[i] + files[i] + '.grp' + haloid[i] + '.haloid.dat',stepnames,stephaloids,format='(A,A)'
            k = where(stepnames EQ files[i] + '.' + steps[valid[j],i] + '/' + files[i] + '.' + steps[valid[j],i])            
            validids[j] = stephaloids[k]
        ENDFOR
        halo_info,$
         dirs[i] + files[i] + '.' + steps[valid,i] + '/' + files[i] + '.' + steps[valid,i],$
          dirs[i] + files[i],$
          haloid[i],$
          validids,$
          massunits[i],$
          distunits[i],$
          debug = debug
    ENDFOR
ENDIF

;------------------------------------- Metallicity ---------------------------------
IF do_metal THEN BEGIN
    uniq = uniq(dirs)
    dirs_uniq = dirs[uniq]
    files_uniq = files[uniq]
    steps_uniq = steps[*,uniq]
    FOR i = n_elements(uniq) - 2, n_elements(uniq) - 1 DO BEGIN
        valid = where(steps_uniq[*,i] NE '0')
        FOR j = 0, n_elements(valid) - 1 DO BEGIN
            cd,dirs_uniq[i] + '/' +  files_uniq[i] + '.' + steps_uniq[valid[j],i]
            metals = mzr(files_uniq[i] + '.' + steps_uniq[valid[j],i]) 
            mwrfits,metals,files_uniq[i] + '.' + steps_uniq[valid[j],i] + '.metals.fits',/create
        ENDFOR
    ENDFOR
ENDIF

IF write_metal THEN BEGIN
    openw,1,'~/Plots/profile/metal.txt'
    FOR i = 0, n_elements(dirs) - 1 DO BEGIN
        valid = where(steps[*,i] NE '0')
        validids = strarr(n_elements(valid))
        FOR j = 0, n_elements(valid) - 1 DO BEGIN
            readcol,dirs[i] + files[i] + '.grp' + haloid[i] + '.haloid.dat',stepnames,stephaloids,format='(A,A)'
            k = where(stepnames EQ files[i] + '.' + steps[valid[j],i] + '/' + files[i] + '.' + steps[valid[j],i])            
            validids[j] = stephaloids[k]
            metal_data = mrdfits(dirs[i] + files[i] +  '.' + steps[valid[j],i] + '/' + files[i] + '.' + steps[valid[j],i] + '.metals.fits',1)
            k = where(metal_data.grp EQ validids[j])
            z_sfr = (2.09*10^(metal_data[k].ox_sfr - 12) + 1.06*metal_data[k].fe_sfr)/zsolar
            printf,1,files[i],' ',haloid[i],' ',steps[valid[j],i],z_sfr
            print,   files[i],' ',haloid[i],' ',steps[valid[j],i],z_sfr
            stop
        ENDFOR
    ENDFOR
    close,1
ENDIF

;--------------------------Global SK law ----------------------------
IF  plot_schmidtlaw_global_obs THEN BEGIN
    angle = 45.0
    intq = 33
    cameras = 15
    center = [0,0]
    j = (size(steps))[1] - 1;7

    IF (KEYWORD_SET(outplot)) THEN BEGIN
        device,filename = outplot + '_SK.eps',/encapsulated,/color,/times,ysize=xsize,xsize=xsize,bits_per_pixel= 8
    ENDIF ELSE BEGIN
        set_plot,'x'
        window,1
    ENDELSE
    IF KEYWORD_SET(color) THEN loadct,39
    xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
    ysigma=2.5e-4*xsigma^1.4
 
    readcol,datafile_prefix + 'HIcubes/ks98.dat',name,D,logHI,logH2,logH,logSFR,tdyn,co,HI,halpha ;,format='(A9D)'
    readcol,datafile_prefix + 'HIcubes/uavpl.dat',logHdarf,logSFRdwarf

    loadct,ctables[0]    
    plot,alog10(xsigma),alog10(ysigma),ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,xrange = [-0.5,2.5], yrange = [-4,-0.5], xtitle=textoidl('Log \Sigma')+"!lH!n [M"+sunsymbol()+" pc!u-2!n]"   
    IF keyword_set(color) THEN loadct,39 ELSE loadct,0
    oplot,logH,logSFR,psym = 2,color = obscolor,symsize = 2
    oplot,logHdarf,logSFRdwarf,psym = 1,color = obscolor,symsize = 2 

    valid = where(steps[j,*] NE '0')
    validids = strarr(n_elements(valid))
    data=fltarr(2,n_elements(valid))
    FOR i = 0, n_elements(valid) - 1 DO BEGIN
        readcol,dirs[valid[i]] + files[valid[i]] + '.grp' + haloid[valid[i]] + '.haloid.dat',stepnames,stephaloids,format='(A,A)'
        k = where(stepnames EQ files[valid[i]] + '.' + steps[j,valid[i]] + '/' + files[valid[i]] + '.' + steps[j,valid[i]])
        validids[i] = stephaloids[k]
        filename = dirs[valid[i]] + files[valid[i]] + '.' + steps[j,valid[i]] + '/' + files[valid[i]] + '.' + steps[j,valid[i]] + '.halo.' + validids[i]
        print,filename
;        data[*,i] = schmidtlaw_global_obs(dirs[valid[i]] + '/steps/' + files[valid[i]] + '.' + steps[j,valid[i]] + '.dir/',files[valid[i]] + '.' + steps[j,valid[i]],dirs[valid[i]] + files[valid[i]] + '.param',/useH2,halo_str = validids[i],angle = 45.0,intq = intq, camera = cameras,center = center,verbose = verbose,/intrinsic);,/Halpha,tipsyfile = filebases[i]
        rtipsy,filename + '.std',h,g,d,s
        units = tipsyunits(dirs[valid[i]] + files[valid[i]] + '.param',/silent)
        read_tipsy_arr,filename + '.HI',h,HI,part = 'gas',type = 'float'
        read_tipsy_arr,filename + '.H2',h,H2,part = 'gas',type = 'float'
        H2 = H2*2.0

        IF file_test(dirs[valid[i]] + files[valid[i]] + '.starlog') THEN BEGIN
            read_tipsy_arr,filename + '.iord',h,siord, part = 'star'
            sl = rstarlog(dirs[valid[i]] + files[valid[i]] + '.starlog',/molecularH)
            match,siord,sl.iorderstar,ind1,ind2
            smassform = sl[ind2].massform
        ENDIF ELSE BEGIN
            read_tipsy_arr,filename + '.iord',h,siord, part = 'star',type = 'long'
            rtipsy,dirs[valid[i]] + files[valid[i]] + '.' + steps[j,valid[i]] + '/' + files[valid[i]] + '.' + steps[j,valid[i]],hall,g,d,s,/justhead
            read_tipsy_arr,dirs[valid[i]] + files[valid[i]] + '.' + steps[j,valid[i]] + '/' + files[valid[i]] + '.' + steps[j,valid[i]] + '.iord',hall,siordall, part = 'star',type = 'long'
            read_tipsy_arr,dirs[valid[i]] + files[valid[i]] + '.' + steps[j,valid[i]] + '/' + files[valid[i]] + '.' + steps[j,valid[i]] + '.massform',hall,smassformall,part = 'star',type = 'float'
            match,siord,siordall,ind1,ind2
            smassform = smassformall[ind2]
        ENDELSE

        g.mass = g.mass*units.massunit
        g.x = g.x*units.lengthunit*h.time
        g.y = g.y*units.lengthunit*h.time
        g.z = g.z*units.lengthunit*h.time
        gr = sqrt(g.x*g.x + g.y*g.y + g.z*g.z)
        s.mass = s.mass*units.massunit
        smassform = smassform*units.massunit
        s.x = s.x*units.lengthunit*h.time
        s.y = s.y*units.lengthunit*h.time
        s.z = s.z*units.lengthunit*h.time
        sr = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)
        currenttime = max(s.tform*units.timeunit)
        deltat = 100*1e6
        recentSF = where(s.tform*units.timeunit GT currenttime - deltat)
        y = weighted_histogram(sr[recentSF],weight = smassform[recentSF],locations = x,nbins = 100,/cum)
        temp = min(abs(y[uniq(y)]/max(y) - 0.8),fiti)
        sfrad = spline(y[(uniq(y))[fiti - 1:fiti + 1]]/max(y),x[(uniq(y))[fiti - 1:fiti + 1]],[0.8],0.01)
;        sigmaH = total((HI[where(gr LE sfrad)] + H2[where(gr LE sfrad)])*g[where(gr LE sfrad)].mass)*units.massunit/(!PI*1000.0*sfrad*1000.0*sfrad)
        sigmaH = total((HI[where(gr LE sfrad)])*g[where(gr LE sfrad)].mass)/(!PI*1000.0*sfrad*1000.0*sfrad)
        sigmaSF = total(smassform[where(sr LE sfrad AND s.tform*units.timeunit GT currenttime - deltat)])/(!PI*sfrad*sfrad)/deltat
        print,sfrad,alog10(sigmaH),alog10(sigmaSF)
        rmax = 30
        nbins = rmax
        sfr_prof = prof_array(sr[recentSF],smassform[recentSF]/deltat,rmin = 0, nbins = nbins, rmax = rmax)
        HI_prof = prof_array(gr,g.mass,weight = HI, rmin = 0, nbins = nbins, rmax = rmax)
        H2_prof = prof_array(gr,g.mass,weight = H2, rmin = 0, nbins = nbins, rmax = rmax)
        print,max(alog10(HI_prof.sd/1e6)),max(alog10(sfr_prof.sd))
 
        data[*,i] = [sigmaH,sigmaSF]
        oplot,[alog10(data[0,i]),alog10(data[0,i])],[alog10(data[1,i]),alog10(data[1,i])],psym = symcat(psym[i]),color = colors[i],symsize = symsizes[i] 
        oplot,alog10(HI_prof.sd/1e6),alog10(sfr_prof.sd),psym = symcat(psym[i]),color = colors[i],symsize = symsizes[i]/3.0
    ENDFOR 
 IF (KEYWORD_SET(outplot)) THEN device,/close ELSE stop
ENDIF

;----------------------- Star Formation History ------------------------------
IF sfh THEN BEGIN
    j = 11
    FOR i = 17, n_elements(files) - 1 DO BEGIN
        filename = dirs[i] + files[i] + '.' + steps[j,i] + '/' + files[i] + '.' + steps[j,i] + '.halo.' + haloid[i]
        print,filename
;        data[*,i] = schmidtlaw_global_obs(dirs[valid[i]] + '/steps/' + files[valid[i]] + '.' + steps[j,valid[i]] + '.dir/',files[valid[i]] + '.' + steps[j,valid[i]],dirs[valid[i]] + files[valid[i]] + '.param',/useH2,halo_str = validids[i],angle = 45.0,intq = intq, camera = cameras,center = center,verbose = verbose,/intrinsic);,/Halpha,tipsyfile = filebases[i]
        units = tipsyunits(dirs[i] + files[i] + '.param',/silent)
        IF i NE 7 AND i NE 8 AND i NE 9 THEN BEGIN
            rtipsy,filename + '.std',h,g,d,s,/justhead
            starlog = rstarlog(dirs[i] + files[i] + '.starlog',/molecularH)
            read_tipsy_arr,filename + '.iord',h,iord,part = 'star',type = 'long'
            match,starlog.iorderstar,iord,ind,ind2
;            stop
            sfr,starlog[ind],meansfr,massunit = units.massunit,timeunit = units.timeunit,starlog  = 1,sarray = sarray,tarray = tarray,mint = 1e-6,maxt = 1.38d10,title = files[i] + '.' + steps[j,i] + '.halo.' + haloid[i]
   ;         stop
        ENDIF ELSE BEGIN
            rtipsy,filename + '.std',h,g,d,s
            stop
            sfr,s,meansfr,massunit = units.massunit,timeunit = units.timeunit,massform = units.istarmass/units.massunit,sarray = sarray,tarray = tarray,mint = 1e-6,maxt = 1.38d10,title = files[i] + '.' + steps[j,i] + '.halo.' + haloid[i]
  ;          stop
        ENDELSE
        writecol,dirs[i] + files[i] + '.halo.' + haloid[i] + '.sfh.txt',tarray,sarray
    ENDFOR
ENDIF
END
