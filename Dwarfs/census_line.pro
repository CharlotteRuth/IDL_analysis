pro census_line_master
mass_str = ['12','10']
mass = [12,10]
imf = ['k']
base = ['/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/','/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/1E5R/']
res_str0 = ['5E1R','1E2R','1E3R','1E4R','1E5R','1E6R']
res_str1 = ['8kpc','2kpc','1kpc','500pc','200pc','50pc']
resolution0 = [50,100,1000,10000,100000,1000000]
resolution1 = [2e-2,1e-2,5e-3,2.5e-3,1e-3,2.5e-4]
ymin=[0.01,0.01]
thickness = [3,3,3] ;[3,1,3,1]
linestyle = [0,2,3] ;[0,0,2,1]
xtitle = ['Number of DM Particles','Softening Length ' + textoidl('[R_{vir}]')]
ytit='Mass Fraction'
loadct,39
!P.CHARSIZE = 1.25
!P.thick = 1.5
!X.Charsize = 1
!Y.Charsize = 1
!X.style = 1
!P.FONT = 0

FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
    set_plot,'x'
    !P.MULTI = [0]
    multiplot,/reset
    set_plot,'ps'
    device,filename = '/astro/net/scratch1/christensen/DwarfResearch/results/census_'+imf[imfct]+'.eps',/color,bits_per_pixel=8,/times

    !P.MULTI = [0,2,2,0,0]
    FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN  
            IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00'
            multiplot
            IF (mct eq 0) THEN BEGIN 
                IF (ploti eq 0) THEN census_line, files, resolution0, mass[mct], ytit = ytit, ymin = ymin[ploti],linestyle=linestyle,thickness=thickness,ytickname=[' ','0.1','0.2','0.3','0.4'] ELSE census_line, files, resolution1, mass[mct], ymin = ymin[ploti],linestyle=linestyle,thickness=thickness 
            ENDIF ELSE BEGIN
                IF (ploti eq 0) THEN  census_line, files, resolution0, mass[mct], xtit = xtitle[ploti], ytit = ytit, ymin = ymin[ploti],linestyle=linestyle,thickness=thickness,ytickname=['0.0','0.1','0.2','0.3','0.4'] ELSE census_line, files, resolution1, mass[mct], xtit = xtitle[ploti], ymin = ymin[ploti],linestyle=linestyle,thickness=thickness,xtickname = [textoidl('10^{-3}'),textoidl('10^{-4}')] 
            ENDELSE
            IF(ploti eq 1 AND mct eq 0) THEN legend,['Stellar','Accreting','Blowaway'],linestyle = linestyle[0:2],thick=thickness[0:2],/right,/bottom,pspacing = 1.5*1.25,charsize = 1.0
            ;legend,['Accreted','Stellar','Blowout','Blowaway'],linestyle = linestyle,thick=thickness,/right,/bottom,pspacing = 1.5*1.25,charsize = 1.0
        ENDFOR
    ENDFOR
    xyouts,[3000,10100,3000,10100],[11100,11100,6000,6000],[textoidl('10^{12}M')+sunsymbol(),textoidl('10^{12}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol(),textoidl('10^{10}M')+sunsymbol()],charsize = 1.0,/device
    device,/close
ENDFOR
cd,'/astro/users/christensen/code/Dwarfs'
END

pro census_line,files,resolutionshort,mass,xtit = xtit,ytit = ytit, title = title,ymin = ymin,linestyle=linestyle,thickness=thickness,_EXTRA = _EXTRA
col=[64,0,254]
m=2.325e5 ;Mass unit
t=1e9 ;time unit
bin=1e8 ;time bin size
resolutions=[resolutionshort[0],resolutionshort,resolutionshort[N_ELEMENTS(resolutionshort)-1]]
xrange = [resolutions[0],resolutions[N_ELEMENTS(resolutions)-1]]
length=N_ELEMENTS(files)
times = ['200','250','300']

starfrac=fltarr(length)
accfrac=fltarr(length)
bofrac=fltarr(length)
awaymassfrac=fltarr(length)
notfrac=fltarr(length)

mtot=fltarr(length)
mtol=fltarr(length)
ms=fltarr(length)

;!p.multi=[0,3,2]
for i=0,length-1 do begin
    print,files[i]
    FOR it = 0,N_ELEMENTS(times) - 1 DO BEGIN
        file = files[i]+ times[it]
    
        mvir=(10.^mass)/m       ;Viral Mass in System Units
        rvir=(3.*mvir*m/(4.*!PI*200.*278*0.7^2))^(1./3.) ;Viral Radius
        
        s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
        rtipsy,file,h,g,d,s
        gr = (g.x^2. + g.y^2. + g.z^2.)^(0.5) ;Gas Radius
        dr = (d.x^2. + d.y^2. + d.z^2.)^(0.5) ;Dark Matter Radius
        sr = (s.x^2. + s.y^2. + s.z^2.)^(0.5) ;Stellar Radius
        gv = (g.vx^2. + g.vy^2. + g.vz^2.)^(0.5) ;Gas Velocity 
        gv2 = (g.vx^2. + g.vy^2. + g.vz^2.) ;Gas Velocity Squared
; Just use z component since that's the way the wind is blowing
        gdot = g.vz*g.z         ;(g.vx*g.x + g.vy*g.y + g.vz*g.z)
        gdens = g.dens*9.4952922e-3 ;Gas Denisty in System Units
        gmr=gr
        be = 0.5*gv2 + g.phi    ;mvir/rvir
        vesc = (2*mvir/rvir)^(0.5) ;Escape Velocity

        accind = where (gdens GT 1e-2 AND g.tempg LT 3e4 AND gdot LE 0, naccind) ;(accretted gas)
        outind = where(gdot GT 0 AND g.zmetal GT 0 AND be LE 0 AND g.vz GT 0.1*vesc, noutind) ;(blow out)
        awayind = where(be GT 0, nawayind, comp=boundind) ;blown away
        
        mstar=total(s.mass)
        mgas=total(g.mass)
        mdark=total(d.mass)
        ms[i]=mstar*m

        starfrac[i] = starfrac[i] + mstar/(mstar+mgas) ;starfrac = fraction of gas as stars
        if (naccind GT 0 ) then accfrac[i] = accfrac[i] + total(g[accind].mass)/(mstar+mgas)$
        else accfrac[i] = accfrac[i]; accfrac = fraction of gas that is falling toward the nucleus  
;        dummy =  total(accind)
;        if dummy eq -1 then accfrac[i] = accfrac[i] $
;        else accfrac[i] = accfrac[i] + total(g[accind].mass)/(mstar+mgas) ; accfrac = fraction of gas that is falling toward the nucleus
        if (noutind GT 0 ) then bofrac[i] = bofrac[i] + total(g[outind].mass)/(mstar+mgas)$
        else bofrac[i] = bofrac[i] ;bofrac = Fraction of gas that will escape perminantly (blow out)
        if (nawayind GT 0 ) then awaymassfrac[i] = awaymassfrac[i] + total(g[awayind].mass)/(mstar+mgas)$
        else awaymassfrac[i] = awaymassfrac[i];nbeind = Fraction of gas beyond mvir/rvir (blow away)        


    ENDFOR 
    starfrac[i] = starfrac[i]/N_ELEMENTS(times)
    accfrac[i] = accfrac[i]/N_ELEMENTS(times)
    bofrac[i] = bofrac[i]/N_ELEMENTS(times)
    awaymassfrac[i] = awaymassfrac[i]/N_ELEMENTS(times)

    print,"Stellar       (starfrac): ",starfrac[i]
    print,"Accreted       (accfrac): ",accfrac[i]
    print,"Blown Out       (bofrac): ",bofrac[i]
    print,"Blown Away (outmassfrac): ",awaymassfrac[i]
    notfrac[i]=1-(starfrac[i]+accfrac[i]+awaymassfrac[i]+bofrac[i])
    mtot[i]=(mstar+mgas+mdark)*m
    mtol[i]=(mstar+mgas+mdark)/mstar

    IF (0) THEN BEGIN
        loadct,39
        v =  SQRT(g.vx*g.vx + g.vy*g.vy + g.vz*g.vz)
        gas_distro = histogram(v/vesc*gdot/ABS(gdot),nbins = 100,locations = vaxes)
        gas_distro_v = histogram(g.vz/vesc,nbins = 100,locations = vaxes_v)
        gas_distro_acc = histogram(v[accind]/vesc*gdot[accind]/ABS(gdot[accind]),nbins = 100,locations = vaxes_acc)
        gas_distro_out = histogram(v[outind]/vesc*gdot[outind]/ABS(gdot[outind]),nbins = 100,locations = vaxes_out)
        IF (nawayind GT 0) THEN gas_distro_away = histogram(v[awayind]/vesc*gdot[awayind]/ABS(gdot[awayind]),nbins = 100,locations = vaxes_away)        


        plot,vaxes,gas_distro,/ylog,yrange = [1e0,1e6],xrange=[-2,2]
        oplot,vaxes_v,gas_distro_v,linestyle = 2
        oplot,vaxes_acc,gas_distro_acc,linestyle = 1,color = 240
        oplot,vaxes_out,gas_distro_out,linestyle = 1,color = 150
        IF (nawayind GT 0) THEN oplot,vaxes_away,gas_distro_away,linestyle = 1,color = 80
        stop
    ENDIF
ENDFOR

;fbind=[0,1,2,4]
;zeroIND=where(outmassfrac GT 0)
;if (zeroIND[0] GT -1) then zero=min(outmassfrac[where(outmassfrac GT 0)]) else zero=min(starfrac[where(starfrac GT 0)])
;if (min(outmassfrac) EQ 0) then outmassfrac[where(outmassfrac EQ 0)] = zero

;ymin=MIN(outmassfrac)
 
plot,resolutionshort,awaymassfrac,xrange=xrange,/xlog,xstyle=1, $
  yrange=[ymin,0.4],ystyle=0, /nodata,xtit=xtit,ytit=ytit, $
  xmargin=[6.5,1.4],ymargin=[3.5,3.5],_EXTRA = _EXTRA
;, $   title=title + ', Galaxy Mass: 10^'+massstr[i_m]+' M'+sunsymbol()

;Stellar Mass fraction (Blue)
oplot,resolutionshort,starfrac, linestyle = linestyle[0],thick=thickness[0]

;Accreted Matter (Purple)
oplot,resolutionshort,accfrac, linestyle = linestyle[1],thick=thickness[1]

;;Blow Out fraction (Orange)
;;oplot,resolutionshort,bofrac,linestyle = linestyle[2],thick=thickness[2]

;Blow Away (Red)
oplot,resolutionshort,awaymassfrac, linestyle = linestyle[2],thick=thickness[2] ;;linestyle = linestyle[3],thick=thickness[3]

; legend,['Accreted Matter','Stellar Mass','Mass Blown Out','Mass Blown Away'],linestyle = [0,0,2,2], color = [0,100,0,100]
; legend,['Accreted Matter','Stellar Mass','Mass Blown Away'],linestyle = [0,0,2], color = [0,100,100]
END
