PRO tdyn_master
files = ['h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK.starlog']
files = ['h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK.starlog',$
         'h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_noH2SF.starlog',$
         'h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_noJeans.starlog',$
         'h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_flatLW.starlog', $
         'h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_lowCstar.starlog', $
         'h516.cosmo25cmb.1536g6MbwK/h516.cosmo25cmb.1536g6MbwK.starlog']

c_star = [0.1,0.1,0.1,0.1,0.01,0.1]
msol_per_sysmass       = 2.310e15
kpc_per_syslength        = 25000.
key = ['Standard','No H2 SF','No Jeans','Flat LW field','Low Cstar','Fabio']

tdyn,files,c_star,msol_per_sysmass,kpc_per_syslength,key
END

PRO tdyn,files,c_star,msol_per_sysmass,kpc_per_syslength,key
close,/all
loadct,39
set_plot,'x'
base = '~/Scratch2/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
set_plot,'ps'

cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
grav_const = 1.10807758d-31 ;in cm^3 amu^-1 s^-2
s_per_year = 3.1557e7
dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
timeunit=SQRT((kpc_per_syslength*3.086d21)^3/(6.67d-8*msol_per_sysmass*1.99d33))/(3600.*24.*365.24)
maxtime = 1.4e+10
colors = (findgen(N_ELEMENTS(files)) + 1)/N_ELEMENTS(files)*254
nel = 500.0
redshift = REVERSE(findgen(nel)/nel*15.0)
lbtime = maxtime - wmap3_lookback(redshift)
nel = 1.0;28.0
tbins = (findgen(nel) + 1.0)/nel*maxtime
binsize = maxtime/nel

FOR it=0,N_ELEMENTS(tbins)-1.0 DO BEGIN
device,/color,bits_per_pixel=8,filename=base+'tdyn_sf.eps',/times,ysize=5,xsize=7,/inch
    !p.multi = [0,2,1]
    close,/all
    FOR i = 0,N_ELEMENTS(files) - 1 DO BEGIN
        star_info = rstarlog(base + files[i])
        star_info.timeform = star_info.timeform*timeunit
        IF ((where(star_info.timeform lt tbins[it]))[0] ne -1) THEN BEGIN
            star_info = star_info[where(star_info.timeform lt tbins[it])] 
            zform = spline(lbtime,redshift,star_info.timeform)
            IF ((where(zform lt 0))[0] ne -1) THEN zform[where(zform lt 0)] = 0
            aform = 1/(zform + 1)
            star_info.rhoform = star_info.rhoform*dens_convert;/aform/aform/aform
            tdynam = 1/SQRT(4.0*!PI*star_info.rhoform*grav_const)/s_per_year
            y = histogram(alog10(tdynam),nbins = 50,locations = x,min = 5.5,max = 8)/DOUBLE(N_ELEMENTS(tdynam))
            IF i eq 0 then plot,x,y,psym = 10,xtitle = 'LOG(Dynamical Time)',yrange = [0,0.2] ;0.12
            oplot,x,y,psym = 10,color = colors[i]
        ENDIF
    ENDFOR
    oplot,[alog10(6.65e6),alog10(6.65e6)],[0,1]
    legend,key,color = colors,linestyle = fltarr(N_ELEMENTS(colors))
    
    FOR i = 0,N_ELEMENTS(files) - 1 DO BEGIN
        star_info = rstarlog(base + files[i])
        star_info.timeform = star_info.timeform*timeunit
        IF ((where(star_info.timeform lt tbins[it]))[0] ne -1) THEN BEGIN
            star_info = star_info[where(star_info.timeform lt tbins[it])]
            y = histogram(star_info.timeform,min = 0,max = maxtime,nbins = nel,locations = x)*MEAN(star_info.massform)*msol_per_sysmass/binsize
;            IF i eq 0 then plot,x/1e10,y,psym = 10,xtitle = 'Time [yr]',ytitle = 'SFR',yrange = [0,0.4],xrange = [0,maxtime/1e10]
;            oplot,x/1e10,y,psym = 10,color = colors[i]
            
        zform = spline(lbtime,redshift,star_info.timeform)
        IF ((where(zform lt 0))[0] ne -1) THEN zform[where(zform lt 0)] = 0
        aform = 1/(zform + 1)
        star_info.rhoform = star_info.rhoform*dens_convert;/aform/aform/aform
        tdynam = 1/SQRT(4.0*!PI*star_info.rhoform*grav_const)/s_per_year
        prob = c_star[i]/tdynam
         y = histogram(alog10(prob),nbins = 50,locations = x,min = -10.0,max = -6.5)/DOUBLE(N_ELEMENTS(tdynam))
        IF i eq 0 then plot,x,y,psym = 10,xtitle = 'LOG(Prob)' ,yrange = [0,0.2],xrange = [-6.5,-10.0]
        oplot,x,y,psym = 10,color = colors[i]
        ENDIF
    ENDFOR
;    stop
    !p.multi = 0
device,/close
ENDFOR
END
