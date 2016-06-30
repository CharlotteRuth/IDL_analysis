PRO plot_profile,filenames,munit,lunit,filebase = filebase,outplot = outplot,debug = debug,starlog = starlog,galtotal = galtotal

loadct,39
formatplot,outplot = outplot
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
H_per_gm = 5.9753790e+23
grav = 6.67e-8
zsolar  =  0.0130215
dtime = 1e8 ;100 MYR

IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 6.25
    ysize = 10
    fileout = outplot
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 500
    ysize = 800
    fileout = filenames
ENDELSE

IF NOT keyword_set(filebase) THEN filebase = filenames
FOR i = 0, n_elements(filenames) - 1 DO BEGIN
    print,filenames[i]
    tunit=sqrt((lunit[i]*3.086d21)^3/(6.67d-8*munit[i]*1.99d33))/(3600.*24.*365.24)
    rtipsy,filenames[i] + '.std',h,g,d,s

;    readarr,filenames[i] + '.HI',h,HI,part = 'gas',/ascii
;    readarr,filenames[i] + '.H2',h,H2,part = 'gas',/ascii
;    readarr,filenames[i] + '.HeI',h,HeI, part = 'gas',/ascii
;    readarr,filenames[i] + '.HeII',h,HeII, part = 'gas',/ascii
;    readarr,filenames[i] + '.FeMassFrac',h,fe, part = 'gas',/ascii
;    readarr,filenames[i] + '.OxMassFrac',h,ox, part = 'gas',/ascii

    read_tipsy_arr,filenames[i] + '.HI',h,HI,part = 'gas'
    read_tipsy_arr,filenames[i] + '.H2',h,H2,part = 'gas'
;    read_tipsy_arr,filenames[i] + '.HeI',h,HeI, part = 'gas'
 ;   read_tipsy_arr,filenames[i] + '.HeII',h,HeII, part = 'gas'
    read_tipsy_arr,filenames[i] + '.FeMassFrac',h,fe, part = 'gas'
    read_tipsy_arr,filenames[i] + '.OxMassFrac',h,ox, part = 'gas'

    cutpos = strsplit(filenames[i],'/')
    slfilename = strmid(filenames[i],0,cutpos[n_elements(cutpos) - 1] - 6) + 'starlog'
    IF file_test(slfilename) THEN BEGIN
;    IF 1 THEN BEGIN
        read_tipsy_arr,filenames[i] + '.iord',h,siord, part = 'star'
        sl = rstarlog(slfilename,/molecularH)
        match,siord,sl.iorderstar,ind1,ind2
        smassform = sl[ind2].massform*munit[i]
    ENDIF ELSE BEGIN
        readarr,filenames[i] + '.iord',h,siord, part = 'star',type = 'long',/ascii
        cutpos = strsplit(filenames[i],'.')
        rtipsy,strmid(filenames[i],0,cutpos[n_elements(cutpos) - 2] - 1),hall,gall,dall,sall,/justhead
        readarr,strmid(filenames[i],0,cutpos[n_elements(cutpos) - 2] - 1) + '.iord',hall,siordall, part = 'star',type = 'long',/ascii
        read_tipsy_arr,strmid(filenames[i],0,cutpos[n_elements(cutpos) - 2] - 1) + '.massform',hall,smassformall, part = 'star'

        match,siord,siordall,ind1,ind2
        smassform = smassformall[ind2]*munit[i]       
    ENDELSE
    close,/all

    H2= H2*2
    z = 1.06*fe + 2.09*ox
    
    vunit = 100.0* 0.73 * (lunit[i] / 1000.0)/2.894405

    d.x = d.x*lunit[i]*h.time
    d.y = d.y*lunit[i]*h.time
    d.z = d.z*lunit[i]*h.time
    drc = sqrt(d.x*d.x + d.y*d.y + d.z*d.z)
    d.mass = d.mass*munit[i]
    epsmin = min(d.eps)*lunit[i]

    g.x = g.x*lunit[i]*h.time
    g.y = g.y*lunit[i]*h.time
    g.z = g.z*lunit[i]*h.time
    g.vx = g.vx*vunit*h.time
    g.vy = g.vy*vunit*h.time
    g.vz = g.vz*vunit*h.time
    gr = sqrt(g.x*g.x + g.y*g.y)
    grc = sqrt(g.x*g.x + g.y*g.y + g.z*g.z)
    g.mass = g.mass*munit[i]
    rhounit=munit[i] * gm_per_msol * H_per_gm/lunit[i]^3/cm_per_kpc^3
    g.dens = g.dens*rhounit/h.time^3

    s.x = s.x*lunit[i]*h.time
    s.y = s.y*lunit[i]*h.time
    s.z = s.z*lunit[i]*h.time
    s.vx = s.vx*vunit*h.time
    s.vy = s.vy*vunit*h.time
    s.vz = s.vz*vunit*h.time
    s.tform = s.tform*tunit
    sr = sqrt(s.x*s.x + s.y*s.y)
    src = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)
    s.mass = s.mass*munit[i]
    now = max(s.tform)
    currentSF = where((now - s.tform) LT dtime)

    sfe = sfe(g,H2 = H2, HI = HI, cstar = 0.1, tempcut = 1e3, denscut = 0.1)
;    He_cold = HeI + HeII
    He = HI;He_cold
    He[where(z LE 0.1)] = (0.236 + 2.1*z[where(z LE 0.1)])/4.0
    He[where(z GT 0.1)] = (-0.446*(z[where(z GT 0.1)] - 0.1)/0.9 + 0.446)/4.0
;    He_cold = He*(HI + H2)/(1 - He - z)
    HHe = (1 - z) ; He + All (molec, atomic, ionized) H
    Hall = HHe - 4*He ;All (molec, atomic, ionized) H
    coldgas = HI + H2 + 4*He + z

    rmax = 100.0
    nbins = 1000
    dr = rmax/nbins
    all_prof = prof_array([drc,grc,src],[d.mass,g.mass,s.mass],rmin = 0, nbins = nbins, rmax = rmax)
    vcirc = sqrt(all_prof.inclosed*gm_per_msol*grav/all_prof.rbins/cm_per_kpc)/1e5
    vcirc[where(all_prof.rbins GT max([drc,grc,src]))] = 0
    s_prof = prof_array(sr,s.mass,rmin = 0, nbins = nbins, rmax = rmax)
    sfr_prof = prof_array(sr[currentSF],smassform[currentSF]/dtime,rmin = 0, nbins = nbins, rmax = rmax)
    HI_prof = prof_array(gr,g.mass,weight = HI, rmin = 0, nbins = nbins, rmax = rmax)
    H2_prof = prof_array(gr,g.mass,weight = H2, rmin = 0, nbins = nbins, rmax = rmax)
    HII_prof = prof_array(gr,g.mass,weight = (Hall - H2 - HI), rmin = 0, nbins = nbins, rmax = rmax)
    Ox_prof = prof_array(gr,g.mass*coldgas,weight = Ox, rmin = 0, nbins = nbins, rmax = rmax)
    Fe_prof = prof_array(gr,g.mass*coldgas,weight = Fe, rmin = 0, nbins = nbins, rmax = rmax)
    Ox_prof_sfe = prof_array(gr,g.mass*sfe,weight = Ox, rmin = 0, nbins = nbins, rmax = rmax)
    Fe_prof_sfe = prof_array(gr,g.mass*sfe,weight = Fe, rmin = 0, nbins = nbins, rmax = rmax)

    s_max = 5.
    g_max = 5.
    ds = 0.05
    dg = 0.05
    s_halfheight = fltarr(nbins)
    g_halfheight = fltarr(nbins)
    s_disp = fltarr(nbins)
    g_disp = fltarr(nbins)

;------------------------------ Calculate PDF of gas density ---------------------    
    mind = alog10(1)
    pdf = weighted_histogram(alog10(g.dens),weight = (g.mass*coldgas), min = mind,max= 5,nbins = 800,locations = pdf_x);min = -3
    pdfterms = gaussfit(pdf_x,pdf,nterms = 3) ;,measure_errors = sqrt(pdf))
    gauss_funct,pdf_x,pdfterms,pdffit
    IF NOT keyword_set(outplot) AND keyword_set(debug) THEN BEGIN
        multiplot,/reset
        window,1
        plot,pdf_x,pdf,xtitle = 'Log Density [amu/cc]'
        oplot,pdf_x,pdffit,color = 100
        print,pdfterms[1],pdfterms[2]
 
        window,3,xsize= 400,ysize = 400
        plot,s.x,s.z,psym = 3,xtitle = 'X [kpc]',ytitle = 'Z [kpc]',xrange = [-20,20],yrange = [-20,20]
        indcold = where(coldgas GT 0.5)
        oplot,g[indcold].x,g[indcold].z,psym = 3,color = 240
       stop
    ENDIF
 ;    diskfit1 = diskfit(s,rmax = 20)

;    fitrange = where(s_prof.rbins LT 20)
;    init = [3.3842650, 1.0872235, 35.089366, 0.18417694, 1.5] 
;    init = [1e8,2.0,1e9,1,1.5]
;    fits = mpfitfun('exp_sersic',s_prof.rbins[fitrange],s_prof.sd[fitrange],sqrt(s_;prof.sd[fitrange]),init)
;    diskfit = abs(fits[0])*exp(-s_prof.rbins[fitrange]/fits[1])
;    diskplot = -2.5*alog10(diskfit)
;    k = 1.9992* fits[4] - 0.3271
;    bulgefit = abs(fits[2])*exp((-k)*((s_prof.rbins[fitrange]/fits[3])^(1/fits[4])-;1))
;    bulgeplot = -2.5*alog10(bulgefit)
;    plot,s_prof.rbins[fitrange],-2.5*alog10(s_prof.sd[fitrange]),yrange = [-15,-25]
;    oplot,s_prof.rbins[fitrange],bulgeplot,linestyle=1,color=254
;    oplot,s_prof.rbins[fitrange],diskplot,linestyle=2,color=254
;    print,"Disk Scale Length: ",strnsignif(fits[1],3),"[kpc]"
;    u_e= strnsignif(-2.5*alog10(sbfactor*fits[2]),3)
;    r_e = strnsignif(fits[3],3)
;    n= strnsignif(fits[4],2)
    
;--------------------------- Calculates the Disk Scale Length and Disk Scale Height -----------------------
    s_diskheight = weighted_histogram(abs(s.z),weight = s.mass        ,min = 0, max = s_max,nbins = s_max/ds, locations = s_locations,/cum)
    s_diskheight = s_diskheight/max(s_diskheight)
    g_diskheight = weighted_histogram(abs(g.z),weight = g.mass*coldgas,min = 0, max = g_max,nbins = g_max/dg, locations = g_locations,/cum)
    g_diskheight = g_diskheight/max(g_diskheight)
    s_locations = s_locations + ds/2.0
    g_locations = g_locations + dg/2.0
    s_uniq = uniq(s_diskheight)
    g_uniq = uniq(g_diskheight)
    IF min(s_diskheight[where(s_diskheight GT 0)]) GT 0.5 OR n_elements(s_uniq) LT 3 THEN s_halfheightall = 0 ELSE $
      s_halfheightall = spline(s_diskheight[s_uniq],s_locations[s_uniq],0.5) ;*max(s_diskheight))
    IF NOT finite(s_halfheightall) Or s_halfheightall LT 0 THEN s_halfheightall = fit_lin([s_diskheight[s_uniq]],[s_locations[s_uniq]],0.5) ;stop
    IF min(g_diskheight[where(g_diskheight GT 0)]) GT 0.5 OR n_elements(g_uniq) LT 3 THEN g_halfheightall = 0 ELSE $
      g_halfheightall = spline(g_diskheight[g_uniq],g_locations[g_uniq],0.5) ;*max(g_diskheight))
    IF NOT finite(g_halfheightall) OR g_halfheightall LT 0 THEN g_halfheightall = fit_lin([g_diskheight[g_uniq]],[g_locations[g_uniq]],0.5) ;stop

    s_disklength = weighted_histogram(sr,weight = s.mass        ,min = 0, max = s_max,nbins = s_max/ds, locations = s_locations,/cum)
    s_disklength = s_disklength/max(s_disklength)
    g_disklength = weighted_histogram(gr,weight = g.mass*coldgas,min = 0, max = g_max,nbins = g_max/dg, locations = g_locations,/cum)
    g_disklength = g_disklength/max(g_disklength)
    s_locations = s_locations + ds/2.0
    g_locations = g_locations + dg/2.0
    s_uniq = uniq(s_disklength)
    g_uniq = uniq(g_disklength)
    IF min(s_disklength[where(s_disklength GT 0)]) GT 0.5 OR n_elements(s_uniq) LT 3 THEN s_halflengthall = 0 ELSE $
      s_halflengthall = spline(s_disklength[s_uniq],s_locations[s_uniq],0.5) ;*max(s_disklength))
    IF NOT finite(s_halflengthall) OR s_halflengthall LT 0 THEN s_halflengthall = fit_lin([0,s_disklength[s_uniq]],[0,s_locations[s_uniq]],0.5) ;stop
    IF min(g_disklength[where(g_disklength GT 0)]) GT 0.5 OR n_elements(g_uniq) LT 3 THEN g_halflengthall = 0 ELSE $
      g_halflengthall = spline(g_disklength[g_uniq],g_locations[g_uniq],0.5) ;*max(g_disklength))
    IF NOT finite(g_halflengthall) OR g_halflengthall LT 0 THEN g_halflengthall = fit_lin([0,g_disklength[g_uniq]],[0,g_locations[g_uniq]],0.5); stop

    print,'Stellar Disk Length: ',strtrim(s_halflengthall,2),' Scale Height: ',strtrim(s_halfheightall,2)
    print,'Gas Disk Length:     ',strtrim(g_halflengthall,2),' Scale Height: ',strtrim(g_halfheightall,2)

     FOR ir = 0, nbins - 1 DO BEGIN
        inds = where(sr GE ir*dr AND sr LT (ir + 1)*dr)
        indg = where(gr GE ir*dr AND gr LT (ir + 1)*dr)
        s_diskheight = weighted_histogram(abs(s[inds].z),weight = s[inds].mass,                     min = 0, max = s_max,nbins = s_max/ds, locations = s_locations,/cum)
        IF max(s_diskheight) NE 0 THEN s_diskheight = s_diskheight/max(s_diskheight)
        g_diskheight = weighted_histogram(abs(g[indg].z),weight = g[indg].mass*(HI[indg] + H2[indg] + 4*He[indg] + z[indg]),min = 0, max = g_max,nbins = g_max/dg, locations = g_locations,/cum)
        IF max(g_diskheight) NE 0 THEN g_diskheight = g_diskheight/max(g_diskheight)
        s_locations = s_locations + ds/2.0
        g_locations = g_locations + dg/2.0
        s_uniq = uniq(s_diskheight)
        g_uniq = uniq(g_diskheight)
;        IF s_diskheight[(epsmin/ds - 1)] GT 0.5 THEN 
;min(s_diskheight[where(s_diskheight GT 0)]) GT 0.5
        IF min(s_diskheight[where(s_diskheight GT 0)]) GT 0.5 OR n_elements(s_uniq) LT 3 OR inds[0] EQ -1 THEN s_halfheight[ir] = 0 ELSE s_halfheight[ir] = spline(s_diskheight[s_uniq],s_locations[s_uniq],[0.5]) ;*max(s_diskheight))
        IF (min(s_diskheight[where(s_diskheight GT 0)]) GT 0.5 OR n_elements(s_uniq) LT 3 OR NOT finite(s_halfheight[ir]) OR s_halfheight[ir] LT 0) AND NOT (s_diskheight[s_uniq] EQ [0] AND n_elements(s_diskheight[s_uniq]) EQ 1) THEN s_halfheight[ir] = fit_lin([0,s_diskheight[s_uniq]],[0,s_locations[s_uniq]],0.5)
;        IF g_diskheight[(epsmin/dg - 1)] GT 0.5 THEN 
        IF s_halfheight[ir] EQ -1 THEN stop

        IF min(g_diskheight[where(g_diskheight GT 0)]) GT 0.5 OR n_elements(g_uniq) LT 3 OR indg[0] EQ -1 THEN g_halfheight[ir] = 0 ELSE g_halfheight[ir] = spline(g_diskheight[g_uniq],g_locations[g_uniq],[0.5]) ;*max(g_diskheight))
        IF (min(g_diskheight[where(g_diskheight GT 0)]) GT 0.5 OR n_elements(g_uniq) LT 3 OR NOT finite(g_halfheight[ir]) OR g_halfheight[ir] LT 0) AND NOT (g_diskheight[g_uniq] EQ [0] AND n_elements(g_diskheight[g_uniq]) EQ 1) THEN g_halfheight[ir] = fit_lin([0,g_diskheight[g_uniq]],[0,g_locations[g_uniq]],0.5)
        IF g_halfheight[ir] EQ -1 THEN stop

        IF 0 THEN BEGIN
;        IF NOT keyword_set(outplot) AND keyword_set(debug) AND (n_elements(g_uniq) GE 2 OR n_elements(s_uniq) GE 2) THEN BEGIN
;        IF s_halfheight[ir] LT 0 OR g_halfheight[ir] LT 0 OR NOT finite(s_halfheight[ir]) OR NOT finite (g_halfheight[ir]) THEN BEGIN
;        IF n_elements(s_uniq) GT 1 OR n_elements(g_uniq) GT 1 THEN BEGIN
;        IF 0 THEN BEGIN
            window,1
            plot,g_locations,g_diskheight,/ylog,yrange = [0.1,1],title = 'Scale Height: ' + strtrim(ir,2) + ' -- ' + strtrim(ir + dr,2),xrange = [0,max(g_locations)]
            oplot,[g_halfheight[ir],g_halfheight[ir]],[0.1,1],linestyle = 2
            oplot,s_locations,s_diskheight,color = 100
            oplot,[s_halfheight[ir],s_halfheight[ir]],[0.1,1],linestyle = 2,color = 100
            oplot,[0,5],[0.5,0.5],linestyle = 1
            print,'Stellar Scale Height: ',s_halfheight[ir]
            print,'Gaseous Scale Height: ',g_halfheight[ir]
            stop
        ENDIF

;        s_disp[ir] = stdev(s[inds].vz)
        s_n_el = n_elements(inds)
        s_weight = s[inds].mass
        s_weight_mean = total(s[indg].vz*s_weight)/total(s_weight)
        s_disp[ir] = sqrt(total(s_weight*(s[inds].vz - s_weight_mean)^2)/((s_n_el - 1)*total(s_weight)/s_n_el))
;        s_disp[ir] = stdev(s[inds].vz)
        IF s_n_el LE 1 THEN s_disp[ir] = 0
        g_n_el = n_elements(indg)
        g_weight = coldgas[indg]*g[indg].mass
        g_weight_mean = total(g[indg].vz*g_weight)/total(g_weight)
        g_disp[ir] = sqrt(total(g_weight*(g[indg].vz - g_weight_mean)^2)/((g_n_el - 1)*total(g_weight)/g_n_el))
;        g_disp[ir] = stdev(g[indg].vz*(HI[indg] + H2[indg]))
        IF g_n_el LE 1 THEN g_disp[ir] = 0
        IF 0 THEN BEGIN
;        IF NOT keyword_set(outplot) AND keyword_set(debug) AND (g_n_el GT 1 OR s_n_el GT 1) THEN BEGIN
            window,1
            gmax = max(weighted_histogram(g[indg].vz,weight = g_weight,nbins = 100,min = -100,max = 100))
            smax = max(weighted_histogram(s[inds].vz,weight = s_weight,nbins = 100,min = -100,max = 100))
            histogramp,g[indg].vz,weight = g_weight,nbins = 100,min = -100,max = 100,title = 'Vel. Disp.: ' + strtrim(ir,2) + ' -- ' + strtrim(ir + dr,2),yrange = [0,max([gmax,smax])]
            oplot,[g_weight_mean,g_weight_mean],[0,1e6],linestyle = 2
            oplot,[g_weight_mean - g_disp[ir],g_weight_mean - g_disp[ir]],[0,1e10],linestyle = 1
            oplot,[g_weight_mean + g_disp[ir],g_weight_mean + g_disp[ir]],[0,1e10],linestyle = 1
            IF s_n_el GT 1 THEN BEGIN
                histogramp,s[inds].vz,weight = s_weight,nbins = 100,min = -100,max = 100,/overplot,color = 100
                oplot,[s_weight_mean,s_weight_mean],[0,1e6],linestyle = 2,color = 100
                oplot,[s_weight_mean - s_disp[ir],s_weight_mean - s_disp[ir]],[0,1e10],linestyle = 1,color = 100
                oplot,[s_weight_mean + s_disp[ir],s_weight_mean + s_disp[ir]],[0,1e10],linestyle = 1,color = 100
            ENDIF
            print,'Stellar Dispersion: ',s_disp[ir]
            print,'Gas Dispersion:     ',g_disp[ir]
            stop
        ENDIF
        
; -------------------  The log  PSF of the gas density in this bin -----------------------
;        pdf = weighted_histogram(alog10(g[indg].dens),weight = (g[indg].mass*(HI[indg] + H2[indg] + 4*He[indg] + z[indg])), min = -2,max= 5,nbins = 800,locations = pdf_x)
;        plot,pdf_x,pdf,xtitle = 'Log Density [amu/cc]'
;        pdfterms = gaussfit(pdf_x,pdf,nterms = 3);,measure_errors = sqrt(pdf))
;        gauss_funct,pdf_x,pdfterms,pdffit
;        oplot,pdf_x,pdffit,color = 100
;        print,10^pdfterms[1],pdfterms[2]
;        stop


        IF 0 THEN BEGIN
            s_diskheight = weighted_histogram(abs(s[inds].z),weight = s[inds].mass        ,min = 0, max = 5,nbins = 50,  locations = locations)
            g_diskheight = weighted_histogram(abs(g[indg].z),weight = g[indg].mass*coldgas,min = 0, max = 5,nbins = 100, locations = locations)
            s_weights = sqrt(s_diskheight)
            g_weights = sqrt(g_diskheight)
            s_indfit = where(s_diskheight GT 1e3)
            g_indfit = where(g_diskheight GT 1e3)
            IF (n_elements(s_indfit) GE 10 AND n_elements(g_indfit) GE 10) THEN BEGIN
                s_a=[max(s_diskheight),0.3]
                s_fit=curvefit(locations[s_indfit],s_diskheight[s_indfit],s_weights[s_indfit],s_a,function_name='expfit')
                g_a=[max(g_diskheight),0.5]
                g_fit=curvefit(locations[g_indfit],g_diskheight[g_indfit],g_weights[g_indfit],g_a,function_name='expfit')
                plot,locations[s_indfit],s_diskheight[s_indfit],/ylog,yrange = [1e3,1e10],title = strtrim(ir,2) + ' -- ' + strtrim(ir + dr,2)
                oplot,locations[s_indfit],s_fit,linestyle = 2
                
                oplot,locations[g_indfit],g_diskheight[g_indfit],color = 100
                oplot,locations[g_indfit],g_fit,linestyle = 2,color = 100
                print,s_a
                print,g_a
            ENDIF
        ENDIF
    ENDFOR

    IF keyword_set(outplot) THEN  device,filename = fileout[i] + '_prof.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  1,yoffset =  .25,/inch ELSE window,0,xsize = xsize, ysize = ysize

;    !p.multi,[1,1,3]
    multiplot,/reset
    multiplot,[1,4],mTitle = filebase[i]
    plot,s_prof.rbins,s_prof.sd/1e6,/ylog,/xlog,ytitle = 'Surface Density [M' + sunsymbol() + textoidl('/pc^2]'),yrange = [1e-4,1e4]
    oplot,HI_prof.rbins,HI_prof.sd/1e6,linestyle = 2
    oplot,H2_prof.rbins,H2_prof.sd/1e6,linestyle = 1
    legend,['Stars','HI',textoidl('H_2')],linestyle = [0,2,1],/right
    multiplot
    
    IF (where(finite(s_halfheight)))[0] EQ -1 THEN s_halfheight = 0*s_prof.rbins
    plot,s_prof.rbins,s_halfheight,/xlog,ytitle = 'Half-Mass Height [kpc]',yrange = [0,5]
    oplot,s_prof.rbins,g_halfheight,linestyle = 2
    legend,['Stars',textoidl('HI + H_2')],linestyle = [0,2],/left
    multiplot

    plot,s_prof.rbins,s_disp,/xlog,ytitle = 'Velocity Dispersion [km/s]',yrange = [0,150]
    oplot,s_prof.rbins,g_disp,linestyle = 2
    legend,['Stars',textoidl('HI + H_2')],linestyle = [0,2],/left
    multiplot

    plot,Ox_prof.rbins,(2.09*Ox_prof.mean + 1.06*Fe_prof.mean),/ylog,/xlog,xrange = [0.5,100],ytitle = 'Log Metallicity',yrange = [1e-5,1e-1]
    oplot,Ox_prof.rbins,Ox_prof.mean,linestyle = 2
    oplot,Fe_prof.rbins,Fe_prof.mean,linestyle = 1
    legend,['Total','Ox','Fe'],linestyle = [0,2,1],/right
    multiplot,/reset

;    ind_num = where(finite(ox_prof_sfe.mean))
;    plot,Ox_prof_sfe[ind_num].rbins,(2.09*Ox_prof_sfe[ind_num].mean + 1.06*Fe_prof_sfe.mean),/ylog,/xlog,xtitle = 'Radius [kpc]',xrange = [0.5,100],ytitle = 'Log Metallicity (Weighted by SF)',yrange = [1e-5,1e-1]
;    oplot,Ox_prof_sfe[ind_num].rbins,Ox_prof_sfe[ind_num].mean,linestyle = 1
;    oplot,Fe_prof_sfe[ind_num].rbins,Fe_prof_sfe[ind_num].mean,linestyle = 2
;    legend,['Total','Ox','Fe'],linestyle = [0,2,1],/right
;    multiplot,/reset
    
    IF keyword_set(outplot) THEN device,/close ELSE BEGIN 
        print,'Filename ','Redshift ','Virial Mass [Msun]','Dark Matter Mass [Msun]','Stellar Mass [Msun]','Half Stellar Mass Radius [kpc]','Half Stellar Mass Height [kpc]',format='(7A40)'
        print,filebase[i],1.0/h.time - 1.0,total([d.mass,s.mass,g.mass]),total(d.mass),total(s.mass),s_halflengthall,s_halfheightall,format='(7A40)'
        print,'Gas Mass [Msun]','HI Mass [Msun]','H2 Mass [Msun]','Half Cold Gas Radius [kpc]','Half Cold Gas Height [kpc]','Log Gas Density Distro Mean ','Log Gas Density Distro STDEV ','Total Gas Metallicity ','Ox ','Fe ','SFR ',format='(11A40)'
        print,total(g.mass),total(g.mass*HI),total(g.mass*H2),g_halflengthall,g_halfheightall,pdfterms[1],pdfterms[2],total(g.mass*sfe*(2.09*Ox + 1.06*Fe))/total(g.mass*sfe),total(g.mass*sfe*ox)/total(g.mass*sfe),total(g.mass*sfe*fe)/total(g.mass*sfe),total(smassform[currentSF])/dtime,format='(10A40)'
        stop
    ENDELSE

;    ind_num = where(finite(sfe))
    openw,1,fileout[i] + '_prof.txt'
    printf,1,'Filename ','Redshift ','Virial Mass [Msun]','Dark Matter Mass [Msun]','Stellar Mass [Msun]','Half Stellar Mass Radius [kpc]','Half Stellar Mass Height [kpc]',format='(7A40)'
    printf,1,filebase[i],1.0/h.time - 1.0,total([d.mass,s.mass,g.mass]),total(d.mass),total(s.mass),s_halflengthall,s_halfheightall,format='(7A40)'
    printf,1,'Gas Mass [Msun]','HI Mass [Msun]','H2 Mass [Msun]','Half Cold Gas Radius [kpc]','Half Cold Gas Height [kpc]','Log Gas Density Distro Mean ','Log Gas Density Distro STDEV ','Total Gas Metallicity ','Ox ','Fe ','SFR ',format='(11A40)'
    printf,1,total(g.mass),total(g.mass*HI),total(g.mass*H2),g_halflengthall,g_halfheightall,pdfterms[1],pdfterms[2],total(g.mass*sfe*(2.09*Ox + 1.06*Fe))/total(g.mass*sfe),total(g.mass*sfe*ox)/total(g.mass*sfe),total(g.mass*sfe*fe)/total(g.mass*sfe),total(smassform[currentSF])/dtime,format='(11A40)'
    printf,1,'Radius [kpc]  ','V_circ [km/s]','Stellar [M_sun kpc^-2]','Stellar Disp. [km/s]','HI [M_sun kpc^-2]','H2 [M_sun kpc^-2]','HII [M_sun kpc^-2]','Cold Gas Disp. [km/s]','Z            ','Z (SF)       ','Ox           ','Ox (SF)       ','Fe           ','Fe (SF)      ','SFR [M_sun kpc^-2',format='(15A23)'
    printf,1,transpose([[HI_prof.rbins],[vcirc],[s_prof.sd],[s_disp],[HI_prof.sd],[H2_prof.sd],[HII_prof.sd],[g_disp],[2.09*Ox_prof.mean + 1.06*Fe_prof.mean],[2.09*Ox_prof_sfe.mean + 1.06*Fe_prof_sfe.mean],[Ox_prof.mean],[Ox_prof_sfe.mean],[Fe_prof.mean],[Fe_prof_sfe.mean],[sfr_prof.sd]]),format='(2f23,e23,f23,3e23,8f23)'
    close,1
    IF keyword_set(debug) THEN stop
ENDFOR

END
