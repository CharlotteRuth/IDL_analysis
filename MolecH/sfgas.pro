PRO sfgas,files,dens_convert = dens_convert,outplot = outplot,color = color,symbols = symbols, thick = thick,key = key,molecularH = molecularH,linestyle = linestyle,_EXTRA = _EXTRA
n = N_ELEMENTS(files)
formatplot,outplot = outplot
IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0
    bgcolor = 255
    device,filename = outplot+'_sfgas.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    window,0,xsize = 712,ysize = 392
ENDELSE

IF NOT keyword_set(dens_convert) THEN dens_convert = 1
IF NOT keyword_set(molecularH) THEN molecularH = fltarr(n)

IF keyword_set(color) THEN BEGIN
    loadct,39
    IF color[0] EQ 1 THEN  color  = (findgen(n) + 1)*240/n ELSE color = color
    IF NOT keyword_set(thick) THEN thick = fltarr(n) + 4
    IF NOT keyword_set(linestyle) THEN linestyle = fltarr(n)
ENDIF ELSE BEGIN
    loadct,0    
    color = (findgen(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(files)) + 5
    IF NOT keyword_set(thick) THEN thick = reverse(findgen(n)*2)
    IF NOT keyword_set(linestyle) THEN linestyle = reverse(findgen(n)*2)
ENDELSE
!p.multi = [0,2,1]

min = 1
max = 4.5
FOR i=0,n -1 DO BEGIN
    data = rstarlog(files[i],molecularH = molecularH[i])
    print,files[i],thick[i],linestyle[i]
    IF i eq 0 THEN histogramp,alog10(data.rhoform*dens_convert),nbins = 100,xtitle = textoidl('log \rho_{Form} [amu/cc]'),ytitle = textoidl('1/N dN/dlog(\rho_{Form} [amu/cc])'),min = min,max = max,/normalize,thick = thick[i],linestyle = linestyle[i],xtickinterval = 1,_EXTRA=_EXTRA;,xrange = [0,4]
    histogramp,alog10(data.rhoform*dens_convert),nbins = 100,min = min,max = max,/normalize,thick = thick[i],linestyle = linestyle[i],color = color[i],/overplot
ENDFOR

min = 1.5
max = 4
FOR i=0,n -1 DO BEGIN
    data = rstarlog(files[i],molecularH = molecularH[i])
    IF i eq 0 THEN histogramp,alog10(data.tempform),nbins = 100,xtitle = textoidl('log T_{Form} [K]'),ytitle = textoidl('1/N dN/dlog(T_{Form} [K])'),min = min,max = max,/normalize,thick = thick[i],linestyle = linestyle[i];,xrange = [0,4]
    histogramp,alog10(data.tempform),nbins = 100,min = min,max = max,/normalize,thick = thick[i],linestyle = linestyle[i],color = color[i],/overplot
ENDFOR

if keyword_set(key) THEN legend,key,color = color,linestyle = linestyle,thick = thick,/right
if (keyword_set(outplot)) then device,/close else stop
!p.multi = 0
END



PRO sfgas_master,vetinari = vetinari

IF keyword_set(vetinari) THEN BEGIN
   prefix = '~/Data/MolecH/Cosmo/' 
   outplot = '~/Figures/h516.cosmo25cmb.paper_'
ENDIF ELSE BEGIN
   prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
   outplot = '~/plots/h516.cosmo25cmb.paper_'   
ENDELSE

files = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.starlog',$
         'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.starlog'] 
pfile = prefix + 'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.param'
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
units = tipsyunits(pfile)
key = ['DnoH2','DH2']
;sfgas,files,dens_convert = units.massunit * gm_per_msol * amu_per_gm/units.lengthunit^3/cm_per_kpc^3,key = key,molecularH = [0,1],outplot = outplot
sfgas,files,dens_convert = units.massunit * gm_per_msol * amu_per_gm/units.lengthunit^3/cm_per_kpc^3,molecularH = [0,1],outplot = '~/Figures/h516.cosmo25cmb.paperColor_',color = [50,240]
END
