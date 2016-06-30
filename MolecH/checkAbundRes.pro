pro checkAbundRes_master
dir = '/astro/store/student-scratch1/christensen/MolecH/12M/'
files = ['Disk_Collapse_1e5/Disk_Collapse_1e5.00100.scalez1.00010','Disk_Collapse_1e6/Disk_Collapse_1e6.00100.scalez1.00010','Disk_Collapse_1e7/Disk_Collapse_1e7.00100.scalez1.00010']
msol_per_sysmass = 2.362e5
kpc_per_syslength = 1.0
outplot = '/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_'

title = ['MWlr','MWmr','MWhr']
checkAbundRes,dir + files,msol_per_sysmass,kpc_per_syslength,title = title,outplot = outplot
end

pro checkAbundRes,files,msol_per_sysmass,kpc_per_syslength,outplot = outplot,title = title,thick = thick,color = color,linestyle = linestyle
n = N_ELEMENTS(files)

formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0
    l_charsize = 0.75
ENDIF ELSE BEGIN
    l_charsize = 1.0
    fgcolor = 255
ENDELSE

IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  color  = (findgen(N_ELEMENTS(broadband)) + 1)*240/N_ELEMENTS(broadband) else color = color
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n)
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n)
ENDIF ELSE BEGIN
    loadct,0    
    color = (findgen(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n)+3
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = reverse(findgen(n)*2)
ENDELSE
linestyle = [1,2,0]



IF NOT KEYWORD_SET(title) THEN title = strind(n)
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
dens_convert =  msol_per_sysmass*gm_per_msol*amu_per_gm/kpc_per_syslength^3/cm_per_kpc^3 ;massunit*amu_per_gm/molec_weight/lengthunit^3 ;Converts to grams/cm^3

!p.multi = 0


IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps'
    device,filename = outplot+'resdens.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 1
ENDIF ELSE BEGIN    
    set_plot,'x'
    window,0,xsize = 712,ysize = 392
ENDELSE

FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s
;    readarr,files[i]+'.H2',h,H2,/ascii,part = 'gas'
;    readarr,files[i]+'.HI',h,HI,/ascii,part = 'gas'
    readarr,files[i]+'.smoothlength',h,hl,/ascii,part = 'gas'

    gamu = g
    gamu.dens = g.dens*dens_convert/h.time/h.time/h.time
    IF (i eq 0) THEN histogramp,alog10(gamu.dens),min = -1,max = 3,nbins = 100,xtitle = textoidl('log \rho [amu/cc]'),linestyle = linestyle[i],/normalize,ytitle = textoidl('1/N dN/d(log \rho [amu/cc])'),yrange = [0,0.0020]
    histogramp,alog10(gamu.dens),min = -1,max = 3,nbins = 100,/overplot,linestyle = linestyle[i],color = color[i],thick = thick[i],/normalize
ENDFOR
legend,title,linestyle = linestyle,thick = thick,color = color,/right
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

stop
IF KEYWORD_SET(outplot) THEN BEGIN
    device,filename = outplot+'resabund.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
    !p.multi = 0
ENDIF ELSE BEGIN
    stop
    window,1,xsize = 712,ysize = 392
    !p.multi = 0
ENDELSE

;!p.multi = [0,n,1]
multiplot,[n,1],mxtitle = textoidl('Density [amu/cc]'),mxTitSize = 1.5,mxTitOffset = 2
FOR i = 0, n - 1 DO BEGIN
    rtipsy,files[i],h,g,d,s
    readarr,files[i]+'.H2',h,H2,/ascii,part = 'gas'
    readarr,files[i]+'.HI',h,HI,/ascii,part = 'gas'
;    readarr,files[i]+'.smoothlength',h,hl,/ascii,part = 'gas'

    gamu = g
    gamu.dens = g.dens*dens_convert/h.time/h.time/h.time
;    hl = hl*cm_per_kpc*kpc_per_syslength*h.time

    IF i eq 0 THEN plot,gamu.dens,2.0*H2/(2.0*H2 + HI),psym = 3,/xlog,xrange = [1e-1,1e4],yrange = [1e-6,1],title = title[i],ytitle = textoidl('H_2/H'),xtickname = [textoidl('10^{-1}'),textoidl('1'),textoidl('10^{1}'),textoidl('10^{2}'),textoidl('10^{3}'),textoidl('10^{4}')],ytickname = ['0.2','0.4','0.6','0.8','1.0'] $;,/ylog $
    ELSE plot,gamu.dens,2.0*H2/(2.0*H2 + HI),psym = 3,/xlog,xrange = [1e-1,1e4],yrange = [1e-6,1],title = title[i],xtickname = [textoidl(' '),textoidl('1'),textoidl('10^{1}'),textoidl('10^{2}'),textoidl('10^{3}'),textoidl('10^{4}')];,/ylog
    print,TOTAL(2.0*TOTAL(H2)/(2.0*TOTAL(H2) + TOTAL(HI)))
    multiplot
ENDFOR

multiplot,/reset


end
