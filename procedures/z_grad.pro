FUNCTION z_grad, filename, pfile = pfile, scale_r = scale_r
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
grav = 6.6725985e-8 ;am^3/g/s
s_per_yr = 3.15569e7

rtipsy,  filename+'.std',h,g,d,s
;readarr, filename+'.FeMassFrac',   h, fe, part = 'gas',/ascii
readarr, filename+'.OxMassFrac',   h, ox, part = 'gas',/ascii
readarr, filename+'.HI',   h, HI, part = 'gas',/ascii
readarr, filename+'.HeI',  h, HeI, part = 'gas',/ascii
readarr, filename+'.HeII', h, HeII, part = 'gas',/ascii
IF file_test(filename+'.H2') THEN readarr, filename+'.H2', h, H2, part = 'gas',/ascii ELSE H2 = HI*0
IF (max(H2) EQ 0) THEN doH2 = 0 ELSE doH2 = 1

a = h.time
IF NOT keyword_set(pfile) THEN BEGIN
    pos = strsplit(filename,'.00',/regex)
    base = strmid(filename,pos[0],pos[1] - pos[0] - 3)
    spawn,'ls ' + base + '*.param',pfile
ENDIF

units = tipsyunits(pfile,scale = a)
gradii = sqrt(g.x*g.x + g.y*g.y + g.z*g.z)*units.lengthunit
;sradii = sqrt(s.x*s.x + s.y*s.y + s.z*s.z)*units.lengthunit

;Weight the gas by the star formation efficiency
IF doH2 THEN BEGIN
    cstar = 0.1
    tempcut = 1e3
    denscut = 0.1
    deltat = 1e6*s_per_yr ;yr
    sfeff = fltarr(h.ngas)
    indsf = where(g.tempg LE tempcut AND g.dens*units.rhounit GE denscut)
    tdyn = 1.0/sqrt(4.0*!PI*grav*g.dens*units.rhounit*gm_per_H)
    IF indsf[0] NE -1 THEN sfeff[indsf] = 1.0 - exp(-1.0*cstar*deltat/tdyn[indsf]*2.0*H2[indsf]/(2.0*H2[indsf] + HI[indsf]))

;    indsf = where(g.tempg LE tempcut)
;    IF indsf[0] NE -1 THEN sfeff[indsf] = 1.0

    sfeff = HI + 2.0*H2
ENDIF ELSE BEGIN
    cstar = 0.1
    tempcut = 1e4
;    denscut = 100.0
    denscut = 10
    deltat = 1e6*s_per_yr ;yr
    sfeff = fltarr(h.ngas)
    indsf = where(g.tempg LE tempcut AND g.dens*units.rhounit GE denscut)
    tdyn = 1.0/sqrt(4.0*!PI*grav*g.dens*units.rhounit*gm_per_H)
    IF indsf[0] NE -1 THEN sfeff[indsf] = 1.0 - exp(-1.0*cstar*deltat/tdyn[indsf])

    sfeff = HI
ENDELSE

IF keyword_set(scale_r) THEN BEGIN
    print,'Scale by R_90'
    pos = strsplit(filename,'.00',/regex)
    base = strmid(filename,pos[0],pos[n_elements(pos) - 1] - pos[0] + 3)
    rscale = enclosedL_R(filename = base,extno = 14,F_num = 4,level = 0.9)
    print,'R_90: ',rscale
    rscale = opticalRadii(filename='h516.cosmo25cmb.3072g14HBWK.00492',B_num = 4)
    print,'R_25: ',rscale
    gradii = gradii/rscale
;    sradii = sradii/rscale

    rmax = 2.5;1.2
    dr = 0.2
ENDIF ELSE BEGIN
    rmax = 8 ;kpc
    dr = 0.1; kpc
ENDELSE

nr = floor(rmax/dr)
grad = replicate({r:0., h_mass:0., h_n:0., ox: 0., fe: 0.},nr)
grad.r = findgen(nr)*dr
;    hrarray  = fltarr(nr)
;    oxrarray = fltarr(nr)
FOR i = 0, nr - 1 DO BEGIN
    ind = where(gradii GE grad[i].r AND gradii LT grad[i].r + dr) 
    print,i,ind[0]
    IF ind[0] NE -1 THEN BEGIN
        subgas = g[ind]
        subgas_mass = g[ind].mass*units.massunit
        subhe_mass = subgas_mass*HeI[ind] + subgas_mass*HeI[ind]
        subhhe_mass = (1 - g[ind].zmetal)*subgas_mass ; He + All (molec, atomic, ionized) H
        subh_mass = subhhe_mass - 4*subhe_mass ;All (molec, atomic, ionized) H
        subh_n =  subh_mass*gm_per_msol/gm_per_H
        subox_n = ox[ind]*subgas_mass*gm_per_msol/2.66d-23
;            subfe_ne = fe[ind]
        subsfeff = sfeff[ind]

        grad[i].h_mass = total(subh_mass)
        grad[i].h_n    = mean(subh_n)
        print,total(subsfeff)
        IF total(subsfeff) NE 0 THEN BEGIN
            grad[i].ox = 12 + alog10(total(subox_n/subh_n*subsfeff)/total(subsfeff))
            grad[i].fe = 0    ;total(subfe_n*subsfeff)/total(subsfeff)
        ENDIF ELSE BEGIN
            grad[i].ox = 0
            grad[i].fe = 0
        ENDELSE
    ENDIF ELSE BEGIN
        grad[i].h_mass = 0
        grad[i].h_n    = 0
        grad[i].ox     = 0
        grad[i].fe     = 0
    ENDELSE
ENDFOR
RETURN,grad
END

;filenames = ['/astro/store/nbody3/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492.halo.1']

;filenames = ['/astro/store/nbody3/christensen/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492.halo.1','/astro/store/nbody3/christensen/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512.halo.1','/astro/store/nbody3/christensen/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.halo.1','/astro/store/nbody3/christensen/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo.1']

;filenames = ['/astro/store/nbody2/abrooks/h277.cosmo50cmb.3072g14HMbwK/h277.cosmo50cmb.3072g14HMbwK.00512/h277.cosmo50cmb.3072g14HMbwK.00512.1']

;base = '/astro/store/nbody3/christensen/'
;filenames = [base + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492.halo.1',base + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512.halo.1',base + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512.halo.1',base + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512.halo.1','/astro/store/nbody2/abrooks/h277.cosmo50cmb.3072g14HMbwK/h277.cosmo50cmb.3072g14HMbwK.00512/h277.cosmo50cmb.3072g14HMbwK.00512.1']
;key = ['h516','h799','h603','h986','h277']
;psym = [4,6,14,16,7]
;z_grad_plot,filenames,key = key, psym = pym,/scale_r

PRO z_grad_plot, filenames, pfiles = pfiles, outfile = outfile, key = key, psym  = psym, color = color, symsize = symsize, ctables = ctables, thicks = thicks, scale_r = scale_r
n = n_elements(filenames)

IF KEYWORD_SET(outfile) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 18
    ysize = 12
    formatplot,/outplot
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 600
    ysize = 400
    formatplot
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if NOT keyword_set(ctables) then ctables = [39,39,39]
    if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(psym) THEN psym = fltarr(n) + 4 ;REVERSE(findgen(n)*2)
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
;   obscolor = fgcolor
;   obssymsize = 1.5
;   obssym = 16
;   obsthick = 2
ENDIF ELSE BEGIN
    loadct,0    
    if NOT keyword_set(ctables) then ctables = [0,0,0]
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(psym) THEN psym = findgen(n) + 4;(findgen(n) + 2)*2   
    IF NOT keyword_set(symsizes) THEN symsizes = fltarr(n) + 2
;   obscolor = 150
;   obssymsize = 1.5
;   obssym = 16
;   obsthick = 2
ENDELSE
IF (KEYWORD_SET(outfile)) THEN BEGIN
   device,filename=outfile+'_z_grad.eps',/color,bits_per_pixel= 8,/times,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2,/encapsul
ENDIF ELSE BEGIN
   window,0,xsize = xsize,ysize = ysize ;392
ENDELSE

FOR i = 0, n - 1 DO BEGIN
    IF keyword_set(pfiles) THEN grad = z_grad(filenames[i],pfile = pfiles[i],scale_r = scale_r) ELSE grad = z_grad(filenames[i],scale_r = scale_r)
    dr = grad[1].r - grad[0].r
    IF i EQ 0 THEN $
      IF keyword_set(scale_r) THEN plot,grad.r + dr/2,grad.ox,xtitle = textoidl('R/R_{90}'), ytitle = '12 + log(O/H)',/nodata,xrange = [0,max(grad.r) + dr],yrange = [7.0,8.2] ELSE plot,grad.r + dr/2,grad.ox,xtitle = 'Radius [kpc]', ytitle = '12 + log(O/H)',/nodata,xrange = [0,max(grad.r) + dr],yrange = [7.6,8.2]
    ind = where(grad.ox NE 0)
    oplot,grad[ind].r + dr/2,grad[ind].ox,psym = -1.0*symcat(psym[i]), color = colors[i], thick = thicks[i], symsize = symsizes[i]
    fit = linfit(grad[ind].r + dr/2,grad[ind].ox,sigma=sigma)
    print,fit,sigma
    xarr = findgen(10)/10.0*10
    oplot,xarr,fit[0] + fit[1]*xarr
ENDFOR
stop
IF keyword_set(outfile) THEN device,/close ELSE stop

END
