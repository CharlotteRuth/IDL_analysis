pro census_line_master
mass_str = ['12','10']
mass = [12,10]

imf = ['ms','k']
;base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/'
;outfile = '/astro/net/scratch1/christensen/DwarfResearch/results/censusNum_'
;directories_num=['5E1R','1E2R', '1E3R','1E4R','1E5R']
;length=N_ELEMENTS(directories_num)
;resolutionshort=[50,100,1000,10000,100000]
;resolutions=[50,resolutionshort,100000]
;xrange=[5e1,1e5]
;xtit='Number of Baryon Particles'
;ytit='Baryonic Mass Fraction'
;title='Mass Distribution at Different Resolutions'

;base = '/astro/net/scratch1/christensen/DwarfResearch/ResSpecTests/'
;outfile = '/astro/net/scratch1/christensen/DwarfResearch/results/censusSpec_'
;directories_spec=['8kpc','2kpc','1kpc','500pc','200pc','50pc']
;resolutionshort=[8000,2000,1000,500,200,50]
;xtit='Spatial Resolution'
;ytit='Baryonic Mass Fraction'
;title='Mass Distribution  (100K particles)'

base = ['../ResMassTests/','../ResSpecTests/1E5R/']
res_str0 = ['5E1R','1E2R','1E3R','1E4R','1E5R']
res_str1 = ['8kpc','2kpc','1kpc','500pc','200pc','50pc']
resulution0 = [50,100,1000,10000,100000]
;resulution1 = [8000,2000,1000,500,200,50]
resolution1 = [2e-2,1e-2,5e-3,2.5e-3,1e-3,2.5e-4]
ymin=[0.01,0.01]
xtit = ['Number of Dark Matter Particles','Force Resolution [pc]']
ytit='Baryonic Mass Fraction'
loadct,39
!P.thick = 1.5
!X.Charsize = 1.25
!Y.Charsize = 1.25
!P.CHARSIZE = 1.25

FOR imfct = 0, N_ELEMENTS(imf) - 1 DO BEGIN
    set_plot,'x'
    !P.MULTI = [0]
;    multiplot,/reset
    set_plot,'ps'
    device,filename = '../results/census_'+imf[imfct]+'.eps',/color,bits_per_pixel=8
    !P.MULTI = [0,2,2,0,0]
    FOR mct = 0, N_ELEMENTS(mass_str) - 1 DO BEGIN
        FOR ploti = 0, N_ELEMENTS(base) -1 DO BEGIN  
            IF (ploti eq 0) THEN files = base[ploti] + res_str0 + '/' + mass_str[mct]+'M_'+imf[imfct]+'/o'+mass_str[mct]+'M_1.00300' ELSE files = base[ploti] + mass_str[mct] + 'M/'+res_str1+'_'+imf[imfct]+'/'+res_str1+'.00300'
;            multiplot
            IF (mct eq 0) THEN BEGIN 
                IF (ploti eq 0) THEN census_line, files, resulution0, mass[mct],xtit = xtit[ploti], ytit = ytit, ymin = ymin[ploti] ELSE census, files, resulution1, mass[mct], xtit = xtit[ploti], ytit = ytit, ymin = ymin[ploti]  
            ENDIF ELSE BEGIN
                IF (ploti eq 0) THEN  census_line, files, resulution0, mass[mct], xtit = xtit[ploti], ytit = ytit, ymin = ymin[ploti] ELSE census, files, resulution1, mass[mct], xtit = xtit[ploti], ytit = ytit, ymin = ymin[ploti]     
            ENDELSE
        ENDFOR
    ENDFOR
    device,/close
ENDFOR
END

pro census_line,files,resolutionshort,mass,xtit = xtit,ytit = ytit, title = title,ymin = ymin
col=[64,0,254]
m=2.325e5 ;Mass unit
t=1e9 ;time unit
bin=1e8 ;time bin size
resolutions=[resolutionshort[0],resolutionshort,resolutionshort[N_ELEMENTS(resolutionshort)-1]]
xrange = [resolutions[0],resolutions[N_ELEMENTS(resolutions)-1]]
length=N_ELEMENTS(files)

outmassfrac=fltarr(length)
bofrac=fltarr(length)
starfrac=fltarr(length)
accfrac=fltarr(length)
notfrac=fltarr(length)
mtot=fltarr(length)
mtol=fltarr(length)
ms=fltarr(length)

;!p.multi=[0,3,2]
for i=0,length-1 do begin
    file = files[i]
    
    mvir=(10.^mass)/m      ;Viral Mass in System Units
    rvir=(3.*mvir*m/(4.*!PI*200.*278*0.7^2))^(1./3.) ;Viral Radius

    s={Mass:0,X:0,Y:0,Z:0,VX:0,VY:0,VZ:0,METALS:0,TFORM:0,EPS:0,PHI:0}
    rtipsy,file,h,g,d,s
    gr = (g.x^2. + g.y^2. + g.z^2.)^(0.5) ;Gas Radius
    dr = (d.x^2. + d.y^2. + d.z^2.)^(0.5) ;Dark Matter Radius
    sr = (s.x^2. + s.y^2. + s.z^2.)^(0.5) ;Stellar Radius
    gv = (g.vx^2. + g.vy^2. + g.vz^2.)^(0.5) ;Gas Velocity 
    gv2 = (g.vx^2. + g.vy^2. + g.vz^2.) ;Gas Velocity Squared
; Just use z component since that's the way the wind is blowing
    gdot = g.vz*g.z             ;(g.vx*g.x + g.vy*g.y + g.vz*g.z)
    gdens = g.dens*9.4952922e-3 ;Gas Denisty in System Units
    gmr=gr
    be = 0.5*gv2 + g.phi        ;mvir/rvir
    vesc = (2*mvir/rvir)^(0.5)  ;Escape Velocity
    beind = where(be GT 0, nbeind, comp=boundind)
    accind = where (gdens GT 1e-2 AND g.tempg LT 3e4 AND gdot LE 0)
    outind = where(gdot GT 0 AND g.zmetal GT 0 AND be LE 0 AND g.vz GT 0.1*vesc, noutind)

    mstar=total(s.mass)
    mgas=total(g.mass)
    mdark=total(d.mass)
    ms[i]=mstar*m
    if (nbeind GT 0 ) then outmassfrac[i]=total(g[beind].mass)/(mstar+mgas)$
    else outmassfrac[i]=0 ;nbeind = Fraction of gas beyond mvir/rvir (blow away)
    
    if (noutind GT 0 ) then bofrac[i]=total(g[outind].mass)/(mstar+mgas)$
    else bofrac[i]=0 ;bofrac = Fraction of gas that will escape perminantly (blow out)
    starfrac[i]=mstar/(mstar+mgas) ;starfrac = fraction of gas as stars
    dummy =  total(accind)
    if dummy eq -1 then accfrac[i]=0 $
    else accfrac[i]=total(g[accind].mass)/(mstar+mgas) ; accfrac = fraction of gas that is falling toward the nucleus

    print,"accfrac: ",accfrac[i]
    print,"starfrac: ",starfrac[i]
    print,"outmassfrac: ",outmassfrac[i]
    notfrac[i]=1-(starfrac[i]+accfrac[i]+outmassfrac[i])
    mtot[i]=(mstar+mgas+mdark)*m
    mtol[i]=(mstar+mgas+mdark)/mstar
endfor

fbind=[0,1,2,4]
zeroIND=where(outmassfrac GT 0)
if (zeroIND[0] GT -1) then zero=min(outmassfrac[where(outmassfrac GT 0)]) else zero=min(starfrac[where(starfrac GT 0)])
if (min(outmassfrac) EQ 0) then outmassfrac[where(outmassfrac EQ 0)] = zero

;ymin=MIN(outmassfrac)
 
plot,resolutionshort,outmassfrac,xrange=xrange,/xlog,xstyle=1, $
  yrange=[ymin,1],/ylog,ystyle=0, /nodata,xtit=xtit,ytit=ytit, $
  xmargin=[6.5,1.4],ymargin=[3.5,3.5]
;, $  title=title + ', Galaxy Mass: 10^'+massstr[i_m]+' M'+sunsymbol()

;Accreted Matter (Purple)
oplot,resolutionshort,accfrac, linestyle = 0,color = 0

;Stellar Mass fraction (Blue)
oplot,resolutionshort,starfrac, linestyle = 0, color = 100

;Blow Out fraction (Orange)
oplot,resolutionshort,bofrac,linestyle = 2, color = 0

;Blow Away (Red)
oplot,resolutionshort,outmassfrac, linestyle =2, color = 100

legend,['Accreted Matter','Stellar Mass','Mass Blown Out','Mass Blown Away'],linestyle = [0,0,2,2], color = [0,100,0,100]

END
