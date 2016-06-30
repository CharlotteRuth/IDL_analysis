pro checkLW,file,dMsolUnit,dKpcUnit,calcfile = calcfile, sigma = sigma,outfile = outfile,center = center,nlw_calc = nlw_calc,changa = changa
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/nobackupp2/crchrist/MolecH/' $
ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/'

gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
H_per_gm = 5.9753790e+23
cmPerKpc = 3.08568025d21
cross_sec = 2.1d-21 ;From Draine and Bertoldi 96.  Same equation as in Gendin 09
zsolar    = 0.0130215
loadct,39

IF NOT keyword_set(center) THEN center = [0,0,0]
formatplot,outplot = outfile
IF keyword_set(outfile) THEN BEGIN
    l_charsize = 0.75
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    l_charsize = 1.0
    fgcolor = 255
    bgcolor = 0
ENDELSE
IF keyword_set(outfile) THEN $
  IF outfile eq '1' THEN outfile = file
IF keyword_set(sigma) THEN outfile = outfile + 'sig'

rtipsy,file,h,g,d,s
readarr,file+".lw",h,lw_all,/ascii
readarr,file+".smoothlength",h,smoothlength,/ascii
smoothlength = smoothlength *dKpcUnit
IF keyword_set(changa) THEN BEGIN
    lw = 10^(lw_all[0:h.ngas - 1] - 2.0*alog10(dKpcUnit) - 2.0*alog10(cmPerKpc)) 
    lw_star = 10d0^(lw_all[n_elements(lw_all) - h.nstar:n_elements(lw_all) - 1])
;    lw_const = 1/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
ENDIF ELSE BEGIN
    lw = lw_all[1:h.ngas]*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
    lw_star = lw_all[n_elements(lw_all) - h.nstar:n_elements(lw_all) - 1]*1d30
;    lw_const = 1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
ENDELSE
;readcol,dir+'MW_disk.00020.lw',lw_mom_all,/silent
;lw_mom = lw_mom_all[1:h.ngas]
;readcol,file+".massform",mf,/silent
;mf2 = mf[where(mf ne 0)]
;mf = mf2[1:n_elements(mf2)-1]

s.x = (s.x - center[0])*dKpcUnit
s.y = (s.y - center[1])*dKpcUnit
s.z = (s.z - center[2])*dKpcUnit
g.x = (g.x - center[0])*dKpcUnit
g.y = (g.y - center[1])*dKpcUnit
g.z = (g.z - center[2])*dKpcUnit
d.x = (d.x - center[0])*dKpcUnit
d.y = (d.y - center[1])*dKpcUnit
d.z = (d.z - center[2])*dKpcUnit

rstar = sqrt(s.x*s.x + s.y*s.y + s.z*s.z);*dKpcUnit
rgas  = sqrt(g.x*g.x + g.y*g.y + g.z*g.z);*dKpcUnit
rdark = sqrt(d.x*d.x + d.y*d.y + d.z*d.z);*dKpcUnit

aveflux = 0
if aveflux then begin
    nel = 20
    maxrad = alog10(100)
    minrad = alog10(0.1)
    radii = 10^(findgen(nel)/nel*(maxrad - minrad) + minrad)
    flux = radii
    flux2 = radii
    flux3 = radii
    nstar = radii
    radii_mass = radii
    x_mass = radii
    y_mass = radii
    z_mass = radii
    for i = 0, nel - 1 do begin
        inds = where(rstar lt radii[i])
        indg = where(rgas lt radii[i])
        indd = where(rdark lt radii[i])    
        radii_mass[i] = 0
        x_mass[i] = 0
        y_mass[i] = 0
        z_mass[i] = 0
        tmass = 0
        if inds[0] ne -1 then begin 
            x_mass[i] = total(s[inds].mass*s[inds].x)
            y_mass[i] = total(s[inds].mass*s[inds].y)
            z_mass[i] = total(s[inds].mass*s[inds].z)
            radii_mass[i] = total(s[inds].mass*rstar[inds])
            tmass = total(s[inds].mass)
        endif
        if indg[0] ne -1 then begin
            x_mass[i] = x_mass[i] + total(g[indg].mass*g[indg].x)
            y_mass[i] = y_mass[i] + total(g[indg].mass*g[indg].y)
            z_mass[i] = z_mass[i] + total(g[indg].mass*g[indg].z)
            radii_mass[i] = radii_mass[i] + total(g[indg].mass* rgas[indg])
            tmass = tmass + total(g[indg].mass)
        endif
        if indd[0] ne -1 then begin 
            x_mass[i] = x_mass[i] + total(d[indd].mass*d[indd].x)
            y_mass[i] = y_mass[i] + total(d[indd].mass*d[indd].y)
            z_mass[i] = z_mass[i] + total(d[indd].mass*d[indd].z)
            radii_mass[i] = radii_mass[i] + total(d[indd].mass*rdark[indd])
            tmass = tmass + total(d[indd].mass)
        endif
        if tmass ne 0 then begin
            x_mass[i] = x_mass[i]/tmass
            y_mass[i] = y_mass[i]/tmass
            z_mass[i] = z_mass[i]/tmass
            radii_mass[i] = sqrt(x_mass[i]^2 + y_mass[i]^2 + z_mass[i]^2)
        endif
        
        if inds[0] ne -1 then begin
            nstar[i] = total(n_elements(s[inds].mass)*mf[inds]) 
            flux[i]  = total(n_elements(s[inds].mass)*mf[inds])/(4.0*3.14159*radii[i]^2) 
            
            dx = s[inds].x - x_mass[i]
            dy = s[inds].y - y_mass[i]
            dz = s[inds].z - z_mass[i]
            d2 = dx*dx + dy*dy + dz*dz
            
            flux2[i] = total(n_elements(s[inds].mass)*mf[inds]/(4.0*3.14159*d2)) 
            flux3[i] = total(n_elements(s[inds].mass)*mf[inds]/(4.0*3.14159*(radii_mass[i])^2))
        endif else begin
            nstar[i] = 0
            flux[i] = 0
            flux2[i] = 0
            flux3[i] = 0
        endelse
    endfor
endif


;--------- Code LW vs radii ----------
IF (keyword_set(outfile)) THEN BEGIN
    device,/close
    device,filename=outfile+'_LW_rad.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
ENDIF ELSE window,0
plot,rgas,lw,/xlog,/ylog,psym = 3,xtitle = 'Radius [Kpc]',ytitle = 'LW';,yrange = [1e-15,100]
stop

;lgas = dbltarr(h.ngas)
;for i =0L, h.nstar - 1 do begin
;    dgas2 = (g.x - s[i].x)^2 + (g.y - s[i].y)^2 + (g.z - s[i].z)^2
;    lgas = lgas + dInitStarMass/(4.0*!PI*dgas2)
;    stop
;endfor
IF NOT keyword_set(nlw_calc) THEN nlw_calc = h.ngas - 1
lgas = dblarr(h.ngas)
print,'Number: ',strtrim(nlw_calc,2)
X_H = 0.7
z_dust = 0.4
g.dens = g.dens*dMsolUnit*gm_per_msol*5.9753790e+23/dKpcUnit^3/cmperkpc^3
averho = 0.1
avemetal = mean(g[where(g.dens gt 0.01)].zmetal)/zsolar

IF NOT keyword_set(calcfile) THEN BEGIN
    IF keyword_set(sigma) THEN tau = averho*avemetal*cross_sec*cmPerKpc*0.7 ELSE tau = 0.0
    FOR i =0L, nlw_calc DO BEGIN
        dgas2 = (g[i].x - s.x)^2 + (g[i].y - s.y)^2 + (g[i].z - s.z)^2
        indd = where(dgas2 LT smoothlength[i]*smoothlength[i])
        IF (indd[0] NE -1) THEN dgas2[indd] = smoothlength[i]*smoothlength[i]
        ind0 = where(dgas2 eq 0,complement = indno0)
        IF (ind0)[0] ne -1 THEN dgas2[ind0] = min(dgas2[indno0])
        lgas[i] = total(lw_star/(4.0*!PI*dgas2*cmPerKpc*cmPerKpc*$
                                 exp(sqrt(dgas2)*tau)))
        IF NOT finite(lgas[i]) then BEGIN
            print,LW_star,dgas2,'NOT FINITE'
            lgas[i] = 0
;            stop
        ENDIF
        IF (i mod 5000) eq 0 THEN print,i 
    endfor
    lgas_all = lw_all 
    lgas_all[1:h.ngas] = lgas/1d30*dKpcUnit*dKpcUnit*cmPerKpc*cmPerKpc
    IF keyword_set(sigma) THEN writecol,file+".lw_calc_sigma",lgas_all ELSE $
      writecol,file+".lw_calc",lgas_all
ENDIF ELSE BEGIN
    read_tipsy_arr,calcfile,h,lgas,part = 'gas'
    lgas = lgas*1d30/dKpcUnit/dKpcUnit/cmPerKpc/cmPerKpc
ENDELSE
oplot,rgas[0:nlw_calc],lgas[0:nlw_calc],psym = 3,color = 60
legend,['LW code','LW calc'],linestyle = [1,1],/bottom,color = [fgcolor,60],charsize = l_charsize


maxh = 1.1 ;kpc
maxr = 8.0 ;kpc
minr = 0.0
diskind = where(abs(g.z) le maxh AND rgas le maxr AND rgas ge minr)
; ----------- Code & Calc LW ------------------------------
IF (keyword_set(outfile)) THEN BEGIN
    device,/close
    device,filename=outfile+'_LW_calc_code.eps.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
ENDIF ELSE window,1
!p.multi = [0,2,1]

; ----------- Code vs Calc LW ------------------------------
;IF (keyword_set(outfile)) THEN BEGIN
;    device,/close
;    device,filename=outfile+'_LW_codeVcalc.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
;ENDIF ELSE window,2
maxlw = 100
minlw = 1e-15
;plot,lw[diskind],lgas[diskind],psym = 3,/ylog,/xlog,xtitle = textoidl('Flux_{LW, Estimate}'),ytitle = textoidl('Flux_{LW}'),xrange = [minlw,maxlw],yrange = [minlw,maxlw]
plot,lw,lgas,psym = 3,/ylog,/xlog,xtitle = textoidl('Flux_{LW, Estimate}'),ytitle = textoidl('Flux_{LW}'),xrange = [minlw,maxlw],yrange = [minlw,maxlw],/nodata
oplot,[minlw,maxlw],[minlw,maxlw]
oplot,lw,lgas,psym = 3,color = 60

; ----------- Code/Calc LW vs radii ------------------------------
;IF (keyword_set(outfile)) THEN BEGIN
;    device,/close
;    device,filename=outfile+'_LW_frac_calc_code.eps.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
;ENDIF ELSE window,1
;plot,rgas[diskind],lw[diskind]/lgas[diskind],psym = 3,/ylog,xrange = [0,maxr],xtitle = 'Radius [Kpc]',ytitle = textoidl('Flux_{LW, Estimate}/Flux_{LW}')
plot,rgas,lw/lgas,psym = 3,/ylog,xrange = [0,maxr],xtitle = 'Radius [Kpc]',ytitle = textoidl('Flux_{LW, Estimate}/Flux_{LW}')
oplot,[minr,maxr],[1,1]
!p.multi = [0]

; ----------- LW vs radius------------------------------
IF keyword_set(ave_radii) THEN BEGIN
    IF (keyword_set(outfile)) THEN BEGIN
        device,/close
        device,filename=outfile+'_LWave_rad.eps',/color,bits_per_pixel= 8,/times
    ENDIF ELSE window,3
    plot,rgas,lw,/xlog,/ylog,psym = 3,xtitle = 'Radius [Kpc]',ytitle = 'LW' ;,yrange = [0.01,1e7],xrange = [0.01,1000],psym = 1
    oplot,rgas[0:nlw_calc],lgas[0:nlw_calc],psym = 3,color = 60
    oplot,radii,flux,linestyle = 1
    oplot,radii,flux2,linestyle = 0
    oplot,radii,flux3,linestyle = 2
    legend,['LW code','LW calc','R = Rmas','R = Rcom - Rstar','R = Rcom'],linestyle = [1,1,1,0,2],/bottom,color = [fgcolor,60,fgcolor,fgcolor,fgcolor],charsize = l_charsize
ENDIF

;-------------- Stellar and gasous profile ---------------
IF 0 THEN BEGIN
    minr = 0.01
    maxr = 10
    ystar = weighted_histogram(rstar,input = s.mass*dMsolUnit, locations = xstar,min = minr,max = maxr,nbins = 200)
    ygas  = weighted_histogram(rgas, input = g.mass*dMsolUnit, locations = xgas,min = minr,max = maxr,nbins = 200)
    IF (keyword_set(outfile)) THEN BEGIN
        device,/close
        device,filename=outfile+'_prof.eps',/color,bits_per_pixel= 8,/times
    ENDIF ELSE window,4
    plot,xstar,ystar,psym = 10,/xlog,/ylog,yrange = [1e6,1e9],xrange = [minr,maxr],xtitle = 'Radius [Kpc]',ytitle = 'Mass'
    oplot,xgas,ygas,psym = 10,linestyle = 1
    legend,['stars','gas'],linestyle = [0,1],/right,charsize = l_charsize  
ENDIF

;-------------- Gal picture, code and calculated ---------------
;print,max(nstar)/(4.0*3.14159*max(radii_mass)^2)
IF (keyword_set(outfile)) THEN BEGIN
    device,/close
    device,filename=outfile+'_LWpic.eps',/color,bits_per_pixel= 8,/times
ENDIF ELSE window,5
!p.multi = [0,2,1]
plot,s.x,s.y,psym = 3,xrange = [-20,20],yrange = [-20,20],title = 'Code';,xrange = [-12,12],yrange = [-12,12]
lwlog = alog10(lw)
lgaslog = alog10(lgas)
minlw = min([lwlog[where(FINITE(lwlog))],lgaslog[where(FINITE(lgaslog))]])
maxlw = max([lwlog[where(FINITE(lwlog))],lgaslog[where(FINITE(lgaslog))]])
color = (lwlog - minlw)/(maxlw - minlw)*200 + 40
for i = 0L, h.ngas - 1 do oplot,[g[i].x,g[i].x],[g[i].y,g[i].y],psym = 3,color = color[i]

plot,s.x,s.y,psym = 3,xrange = [-20,20],yrange = [-20,20],title = 'Calculation';,xrange = [-12,12],yrange = [-12,12]
lgaslog = alog10(lgas)
color = (lgaslog - minlw)/(maxlw - minlw)*200 + 40
for i = 0L, h.ngas - 1 do oplot,[g[i].x,g[i].x],[g[i].y,g[i].y],psym = 3,color = color[i]
!p.multi = 0

IF (keyword_set(outfile)) THEN BEGIN
    device,/close
    device,filename=outfile+'_LWhist.eps',/color,bits_per_pixel= 8,/times
ENDIF ELSE window,6
histogramp,alog10(lw[where(lgas ne 0)]/lgas[where(lgas ne 0)]),xtitle = textoidl('Log (Code Flux_{LW} / Calculated Flux_{LW})'),xrange = [-4,2],nbins = 100,/normalize;,ytitle = textoidl("1/N dN/Log (Code Flux_{LW} / Calc. Flux_{LW})")

IF (keyword_set(outfile)) THEN BEGIN
    device,/close
    device,filename=outfile+'_LWhistlog.eps',/color,bits_per_pixel= 8,/times
ENDIF ELSE window,6

histogramp,alog10(lw[where(lgas ne 0)]/lgas[where(lgas ne 0)]),xtitle = textoidl('Log (Code Flux_{LW} / Calculated Flux_{LW})'),ytitle = textoidl('1/N dN/Log (Code Flux_{LW} / Calc. Flux_{LW})'),xrange = [-4,2],nbins = 100,/normalize,/ylog,yrange = [1e-6,2e-1]
IF (keyword_set(outfile)) THEN BEGIN
    device,/close
ENDIF ELSE stop
end

pro doublecheck
if 0 then begin
    window,2

    y = HISTOGRAM(alog10(lw),binsize = 0.001,locations = x)
    plot,x,y,psym = 10,xrange=[-3,0]


    sub0_303 = where(alog10(lw) le -0.885 AND alog10(lw) ge -0.8878 AND g.x gt 2e-5)
;    sub0_303 = where(alog10(lw) le 0.3049 AND alog10(lw) ge 0.300)


    help,sub0_303
    Vol = (max(g[sub0_303].x) - min(g[sub0_303].x))* $
          (max(g[sub0_303].y) - min(g[sub0_303].y))* $
          (max(g[sub0_303].z) - min(g[sub0_303].z))
 
    starsind = where(s.x gt min(g[sub0_303].x) AND s.x lt max(g[sub0_303].x) AND $
                     s.y gt min(g[sub0_303].y) AND s.y lt max(g[sub0_303].y) AND $
                     s.z gt min(g[sub0_303].z) AND s.z lt max(g[sub0_303].z)) 
    stellarmass = total(mf[starsind])

   print,alog10(stellarmass/(4.0*!PI*vol^(2.0/3.0)))
loadct,39
plot,g[sub0_303].x,g[sub0_303].y,psym = 3
oplot,g[sub0_303].x,g[sub0_303].y,psym = 3,color = 240
;oplot,s.x,s.y,psym = 3
oplot,s[starsind].x,s[starsind].y,psym = 3,color = 50
endif
end


pro checkLW_master
dir = '/astro/net/scratch2/christensen/MolecH/11M/Cloud_1e3/'
file = 'o11M_1.onestar.00001'
file = 'o11M_1.00005'
dMsolUnit       = 2.362e5
dKpcUnit        = 1
;dInitStarMass   = 17.2043010752688

dir = '/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_zsol/rad.psource/'
file = 'MW_disk.00010'
dMsolUnit       = 1.36e17
dKpcUnit        = 1e5
;dInitStarMass   = 2.535817803095e-15   ;3.34632e-14

dir = '/astro/net/scratch2/christensen/MolecH/11M/Cloud_1e5/'
file = 'o11M_00300.00001'
dMsolUnit       = 2.362e5
dKpcUnit        = 1d0

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h277.cosmo50cmb.1536g1MBWKBH/'
file = 'h277.cosmo50cmb.1536g1MBWKBH.00001'
dMsolUnit       = 1.84793e16   
dKpcUnit        = 50000.

dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g7HBWK_00480/steps.JeansSF/h516.cosmo25cmb.1536g8HBWK.JeansSF.00408.00032.dir/'
file = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.00408.00032.halo.1'
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g8HBWK/steps/h516.cosmo25cmb.1536g9HBWK.00408.dir/'
file = 'h516.cosmo25cmb.1536g9HBWK.00408.halo.1'
dMsolUnit       = 2.310e15
dKpcUnit        = 25000.

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/'
file = 'h516.cosmo25cmb.2304g14HBWK.00512.halo.1'
dMsolUnit       = 2.310e15
dKpcUnit        = 25000.

dir = '/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_1e6_comp/Disk_Collapse_1e6_sol_H2/'
file = 'Disk_Collapse.H2.solmet.ch0.00100.00020'
dMsolUnit       = 2.362e5 ;2.310e15
dKpcUnit        = 1; 25000.

dir = '/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_1e6/old_code/'
file = 'Disk_Collapse_1e6.00100.scalez1.00010'
calcfile = dir + file + '.lw_calc'
dMsolUnit       = 2.362e5 ;2.310e15
dKpcUnit        = 1; 25000.

dir = '/home/christensen/Storage1/UW/MolecH/12M/Disk_Collapse_1e5/Disk_Collapse_1e5.00100.LW/'
file = 'Disk_Collapse_1e5.00100.00001'
calcfile = dir + file + '.lw_calc'
dMsolUnit = 2.362e5 
dKpcUnit = 1

cx = 2.159512
cy = 1.099671
cz =-2.671780

;'o11M_1.00005'
cx = 4.991615e-01 
cy =-9.257880e-01 
cz =-2.579842e+00

;'MW_disk.00010'
cx = -5.640686e-07 
cy = -3.050460e-05 
cz = -6.061303e-06

;o11M_00300.00001
cx = -4.692374e-01 
cy = 3.408022e-01 
cz = -6.159824e-01

;h277
cx = 3.719014e-06 
cy = -4.232391e-05 
cz = -9.820476e-06

;Disk_Collpase_1e5.00100.00001
cx = -2.178e-1
cy = 6.3e-1
cz = -9.758665e-2

;checkLW,dir+file,dMsolUnit,dKpcUnit,calcfile = calcfile, sigma=sigma,outfile = outfile,center = [cx,cy,cz],nlw_calc = nlw_calc
outfile = '/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse'
outfile = dir
calcfile = 0
checkLW,dir+file,dMsolUnit,dKpcUnit,calcfile = calcfile,outfile = outfile,center = [cx,cy,cz];,nlw_calc = 5000
end
