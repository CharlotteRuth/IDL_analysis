pro jeans_master
loadct,39
set_plot,'x'
;set_plot,'ps'

colors=[240,210,190,150,100,70]
base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/'
files12 = base + ['1E6R/12M_k/o12M_1.00300','1E5R/12M_k/o12M_1.00300','1E4R/12M_k/o12M_1.00300','1E3R/12M_k/o12M_1.00300','1E2R/12M_k/o12M_1.00300','5E1R/12M_k/o12M_1.00300']
files10 = base + ['1E6R/10M_k/o10M_1.00300','1E5R/10M_k/o10M_1.00300','1E4R/10M_k/o10M_1.00300','1E3R/10M_k/o10M_1.00300','1E2R/10M_k/o10M_1.00300','5E1R/10M_k/o10M_1.00300']

print,'**********************************************************************'
print,'**********************************************************************'
window,0
;device,filename=base +
;'../results/jeans10.eps',/color,bits_per_pixel=8,/times
;[-4,-0.3]
plot,[alog10(32),alog10(32)],[-6,-0.3],xrange = [-2,10],yrange = [-4,-0.3],xtitle = 'Log(Jeans Mass/Praticle Mass)',ytitle = 'Log(Fraction of Gas Particles)',title = '1e10 Solar Mass Halo',ystyle = 1
FOR ct = 0, N_ELEMENTS(files10) - 1 DO jeans,files10[ct],iteration = ct,min = -4
print,'**********************************************************************'
FOR ct = 0, N_ELEMENTS(files10) - 1 DO jeans,files10[ct],iteration = ct,min = -4,/sfgas
legend, ['1e6','1e5','1e4','1000','100','50'],color = colors,linestyle = [0,0,0,0,0,0] 
;device,/close

;window,1
;device,filename=base + '../results/jeans10.eps',/color,bits_per_pixel=8,/times
;plot,[alog10(32),alog10(32)],[-4,-0.3],xtitle = 'Log(Jeans Mass)',ytitle = 'Log(Fraction of Gas Particles)',title = '1e10 Solar Mass Halo',ystyle = 1,min = 5.5,xrange = [5.5,12.5],yrange = [-4,-0.3]
;FOR ct = 0, N_ELEMENTS(files10) - 1 DO jeans,files10[ct],iteration = ct,min = -4,/absolute
;legend, ['1e6','1e5','1e4','1000','100','50'],color = colors,linestyle = [0,0,0,0,0,0] 
;device,/close
;
;stop
print,'**********************************************************************'
print,'**********************************************************************'
window,2
;device,filename=base +
;'../results/jeans12.eps',/color,bits_per_pixel=8,/times
;[-3,-0.5]
plot,[alog10(32),alog10(32)],[-6,-0.5],xrange = [-2,10],yrange = [-3,-0.5],xtitle = 'Log(Jeans Mass/Praticle Mass)',ytitle = 'Log(Fraction of Gas Particles)',title = '1e12 Solar Mass Halo'
FOR ct = 0, N_ELEMENTS(files10) - 1 DO jeans,files12[ct],iteration = ct,min = -3
print,'**********************************************************************'
FOR ct = 0, N_ELEMENTS(files10) - 1 DO jeans,files12[ct],iteration = ct,min = -3,/sfgas
legend, ['1e6','1e5','1e4','1000','100','50'],color = colors,linestyle = [0,0,0,0,0,0] 
;device,/close

;window,3
;device,filename=base + '../results/jeans12.eps',/color,bits_per_pixel=8,/times
;plot,[alog10(32),alog10(32)],[-4,-0.3],xtitle = 'Log(Jeans Mass)',ytitle = 'Log(Fraction of Gas Particles)',title = '1e12 Solar Mass Halo',min = 2.5,yrange = [-3,-0.5],xrange = [2.5,13.5]
;FOR ct = 0, N_ELEMENTS(files10) - 1 DO jeans,files12[ct],iteration = ct,min = -3,/absolute
;legend, ['1e6','1e5','1e4','1000','100','50'],color = colors,linestyle = [0,0,0,0,0,0] 
;device,/close

END

;jeans,'/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g3HBWK_00504/steps_column/h603.cosmo50cmb.2304g3HBWK.00504.00020.dir/h603.cosmo50cmb.2304g3HBWK.00504.00020.halo.1',massunit = 1.84793e16, kpc_per_syslength = 50000.
pro jeans,file,massunit = massunit, kpc_per_syslength = kpc_per_syslength,iteration =iteration,min = min,absolute=absolute,sfgas = sfgas
;jeans,'/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/1E5R/10M_k/o10M_1.00300'
;jeans,'/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/1E6R/10M_k/o10M_1.00300'
colors=[240,210,190,150,100,70]
loadct,39
IF NOT KEYWORD_SET(kpc_per_syslength) THEN kpc_per_syslength = 1.0
IF NOT KEYWORD_SET(massunit) THEN  massunit = 2.362d5 ;system mass unit in solar masses
IF NOT KEYWORD_SET(min) THEN min = -10
IF NOT KEYWORD_SET(iteration) THEN iteration = 0
k = 1.380658e-16 ;Boltzmann constant
mh = 1.673534e-24 ;hydrogen mass in grams
mu = 2.2 ;molecular weight in cloud with 90% H_2 and 10% He
grav = 6.67259e-8 ;Graviational constant in cgs
pi = 3.14159

cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790d23
molec_weight = (0.76*1 + 0.24*4.0)
dens_convert =  massunit * gm_per_msol/kpc_per_syslength^3/cm_per_kpc^3 ;converts density in sysmass/syslength^3 to g/cm^3

xmax = 175
xmin = xmax*(-1.0)
ymax = 175
ymin = ymax*(-1.0)
binsize = 5
nbins = 2*ymax/binsize

;set_plot,'ps', /copy, /interpolate
;device,filename=mstr+'density.eps', /color, bits_per_pixel=8
print,''
print,file
print,'-------------------------'
rtipsy,file,h,g,d,s

j = {X:DOUBLE(0),Y:DOUBLE(0),Z:DOUBLE(0),DENS:DOUBLE(0),TEMPG:DOUBLE(0),MASS:DOUBLE(0),JMASS:DOUBLE(0)}
ind = where(g.tempg lt 15000 AND g.dens*dens_convert gt 0.1)
window,0
plot,[alog10(32),alog10(32)],[-6,-0.3],xrange = [-2,10],yrange = [-4,-0.3],xtitle = 'Log(Jeans Mass/Praticle Mass)',ytitle = 'Log(Fraction of Gas Particles)',title = 'h603.cosmo50cmb.2304g3HBWK',ystyle = 1
IF (NOT KEYWORD_SET(sfgas) OR ind[0] ne -1) THEN BEGIN
    total_particles = N_ELEMENTS(g)
    IF (KEYWORD_SET(sfgas)) THEN BEGIN
        gall = g
        g = g[ind]
    ENDIF
    jeans = REPLICATE(j,N_ELEMENTS(g.x))
    jeans.x = DOUBLE(g.x)
    jeans.y = DOUBLE(g.y)
    jeans.z = DOUBLE(g.z)
    jeans.dens = DOUBLE(g.dens) * dens_convert ;massunit/(lengthunit)^3 in grams per cm^3
    jeans.tempg = DOUBLE(g.tempg)
    jeans.mass = DOUBLE(g.mass*massunit)
    ;jeans.jmass = (5.0*k*g.tempg/(grav*mu*mh))^1.5*(3.0/(4.0*pi*jeans.dens))^0.5 
    ;Jeans mass equation in terms of solar masses
    jeans.jmass = 1.19402e-10*(g.tempg/mu)^1.5*jeans.dens^(-0.5) 
    ;Jeans formula but with collapsed constants
    temp_check = 15000
    dens_check = 0.1*mu*mh
    print,'Jeans Mass (check): ',STRTRIM(LONG64(1.19402e-10*(temp_check/mu)^1.5*dens_check^(-0.5)),2),' Solar Masses'

    ratio = jeans.jmass/jeans.mass
    IF (keyword_set(absolute)) THEN ratio = jeans.jmass
    db = (MAX(alog10(ratio)) - MIN(alog10(ratio)))/50
    ratio_y = alog10(float(histogram(alog10(ratio),min = MIN(alog10(ratio)),binsize = db))/total_particles)
    ratio_x = findgen(N_ELEMENTS(ratio_y))*db+min(alog10(ratio))
    IF ((where(ratio_y lt min))[0] ne -1) THEN ratio_y[where(ratio_y lt min)] = alog10(0)
    if(KEYWORD_SET(sfgas)) then oplot,ratio_x,ratio_y,color = colors[iteration] else oplot,ratio_x,ratio_y,color = colors[iteration], linestyle = 0
    ;print,'Particle Mass:            ',STRTRIM(LONG64(MAX(jeans.mass)),2),' -- ',STRTRIM(LONG64(MIN(jeans.mass)),2),' Solar Masses'
    ;print,'Jeans Mass:               ',STRTRIM(LONG64(MAX(jeans.jmass)),2),' -- ',STRTRIM(LONG64(MIN(jeans.jmass)),2),' Solar Masses'
    ;print,'Jeans Mass/Particle Mass: ',STRTRIM(MAX(ratio),2),' -- ',STRTRIM(MIN(ratio),2)
    print,'Unresolved Jeans Mass: ',FLOAT(N_ELEMENTS(where(ratio lt 32)))/N_ELEMENTS(ratio),' (',N_ELEMENTS(where(ratio lt 32.0)),', ',N_ELEMENTS(ratio),')'
    ind32 = where(ratio lt 32)
    ind10 = where(ratio lt 10)
    ind1 = where(ratio lt 1)

    window,1
    plot,jeans.dens/mh,jeans.tempg,psym = 3,xtitle = 'Density [amu cm^-3]',ytitle = 'Temperature [K]',/xlog,/ylog,yrange = [100,1e7],xrange = [1e-6,1e3]
    oplot,jeans[ind32].dens/mh,jeans[ind32].tempg,psym = 3,color = 200
    oplot,jeans[ind10].dens/mh,jeans[ind10].tempg,psym = 3,color = 220
    oplot,jeans[ind1].dens/mh,jeans[ind1].tempg,psym = 3,color = 240
    stop
    
;image = fltarr(nbins,nbins)
;FOR iy = 0, nbins-2 DO BEGIN
;    IF(iy MOD 10 EQ 0) THEN print,iy
;    FOR ix = 0, nbins-2 DO BEGIN
;        ind = WHERE(((ymin + iy*binsize) LT jeans.y AND (ymin + (iy + 1)*binsize) GT jeans.y) AND ((xmin + ix*binsize) LT jeans.x AND (xmin + (ix + 1)*binsize) GT jeans.x))
;        IF (ind[0] NE -1) THEN image[ix,iy] = MIN(jeans[ind].jmass) ELSE image[ix,iy] = 1e13
;    ENDFOR
;ENDFOR
;image=image[0:68,0:68]
;image=ALOG10(image)
;contour,image,nlevels=7,yrange=[0,68],xrange=[0,68],XSTYLE = 5, YSTYLE = 5,title = 'Jeans Mass',xtitle = 'X (kpc)',ytitle = 'Y (kpc)'

;ind = WHERE(jeans.jmass LE 1e6)
;IF (ind[0] NE -1) THEN oplot,jeans[ind].x/70.0+35.0,jeans[ind].y/70.0+35.0,psym = 2,color = 240

;image_dens = fltarr(nbins,nbins)
;FOR iy = 0, nbins-2 DO BEGIN
;    IF(iy MOD 10 EQ 0) THEN print,iy
;    FOR ix = 0, nbins-2 DO BEGIN
;        ind = WHERE(((ymin + iy*binsize) LT jeans.y AND (ymin + (iy + 1)*binsize) GT jeans.y) AND ((xmin + ix*binsize) LT jeans.x AND (xmin + (ix + 1)*binsize) GT jeans.x))
;        IF (ind[0] NE -1) THEN image_dens[ix,iy] = MEAN(jeans[ind].dens) ELSE image_dens[ix,iy] = 1.59902e-26
;    ENDFOR
;ENDFOR
;image_dens=image_dens[0:68,0:68]
;image_dens=ALOG10(image_dens)
;contour,image_dens,nlevels=7,yrange=[0,68],xrange=[0,68],XSTYLE = 5, YSTYLE = 5,title = 'Gas Density',xtitle = 'X (kpc)',ytitle = 'Y (kpc)',/fill
ENDIF ELSE print,'No Gas Particles'

END
