PRO densRadius,files,dKpcUnit,dMsolunit,outplot = outplot,maxdistance = maxdistance,color = color,symbols = symbols, halos_str = halos_str,verbose = verbose,label = label,ctables = ctables
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
convert2msol = 7.9493857e-21 

n = N_ELEMENTS(files)
mid = CEIL(n/2.0) - 1
;!p.multi = [0,1,n]
if not keyword_set(maxdistance) then maxdistance = 12
;IF NOT KEYWORD_SET(xrange) THEN 
xrange = [0,maxdistance]
IF KEYWORD_SET(outplot) THEN BEGIN
    bgcolor = 255
    fgcolor = 0 
    obscolor = 120
ENDIF ELSE BEGIN
    bgcolor = 0
    fgcolor = 255
    obscolor = 200
ENDELSE
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = color
;    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
;    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) ;REVERSE(findgen(n)*2)
    IF NOT KEYWORD_SET(symbols) THEN symbols = fltarr(n) + 16
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n) + 39
    obscolor = fgcolor
ENDIF ELSE BEGIN
    loadct,0    
    colors = (findgen(n) + 1)*fgcolor;  fltarr(n) + 5
;    IF NOT KEYWORD_SET(thicks) THEN thicks = (findgen(n) + 1)*6/n - 1
;    IF NOT KEYWORD_SET(linestyle) THEN linestyle = REVERSE(findgen(n)*2)   
    IF NOT KEYWORD_SET(symbols) THEN symbols = fltarr(n) + 16 ;(findgen(n)+2)*2
    IF NOT keyword_set(ctables) THEN ctables = fltarr(n)
ENDELSE
IF NOT KEYWORD_SET(halos_str) THEN halos_str = strarr(n) + '1' 
IF KEYWORD_SET(outplot) THEN device,filename=outplot + '_densR.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 ELSE window,0

nbins1 = 20.0
bind1 = maxdistance/nbins1*cm_per_kpc
radiibins = (findgen(nbins1)+1)*bind1
sd = fltarr(nbins1,n)
r25 = fltarr(n)

!X.MARGIN = [12,12]
!Y.MARGIN = [6,4]
multiplot,[1,n],/doyaxis,mtitle = label;,mxtitle = 'Radius [kpc]';,mytitle = 'Log Density [amu/cc]',mxtitle = 'Radius [kpc]';,mytitle = 'Log Density [amu/cc]';,/doxaxis,/doyaxis;gap=0.04,
FOR i = 0, n - 1 DO BEGIN
    r25[i] = opticalRadii(filename = files[i],halo_str = halos_str[i])
    rtipsy,files[i]+ '.halo.' + halos_str[i],h,g,d,s
    read_tipsy_arr,files[i]+ '.halo.' + halos_str[i] + '.HI',h,HIfrac,part = 'gas';,/ascii
    IF (FILE_TEST(files[i]+'.halo.' + halos_str[i] + '.H2"')) THEN BEGIN
        read_tipsy_arr,files[i]+ '.halo.' + halos_str[i] + '.H2',h,H2frac,part = 'gas';,/ascii
        frac = (HIfrac + 2.0*H2frac) ;1.36 incorporates He
    ENDIF ELSE frac = HIfrac*1.36
    gmass = g.mass*frac
    gdens = g.dens*frac
    dens_convert =  dMsolunit[i] * gm_per_msol * 5.9753790e+23/dKpcUnit[i]^3/cm_per_kpc^3/h.time^3
    radius = SQRT(g.x*g.x + g.y*g.y + g.z*g.z)*dKpcUnit[i]*h.time
    disk = where(g.tempg lt 5e4)

    FOR j = 0, nbins1-1 DO BEGIN
        ind = where(radius*cm_per_kpc LE radiibins[j] AND radius*cm_per_kpc GE radiibins[j] - bind1)
        radmid = radiibins[j] + bind1/2
        area = !PI*(radiibins[j]^2 - (radiibins[j] - bind1)^2)
        mass = TOTAL(gmass[ind])*dMsolunit[i] * gm_per_msol * 5.9753790e+23
        sd[j,i] = mass/area
    ENDFOR

    loadct,0
    if i eq n -1 THEN plot,radius[disk],alog10(g[disk].dens*dens_convert*frac),yrange = [-1.99,3],xrange = xrange,psym = 3,xtitle = 'Radius [kpc]',YSTYLE=8,/nodata,xstyle = 1 $ ;,mytitle = 'Density [amu/cc]'
    ELSE IF i eq mid THEN plot,radius[disk],alog10(g[disk].dens*dens_convert*frac),yrange = [-1.99,3],xrange = xrange,psym = 3,YSTYLE=8,ytitle = textoidl('log \rho [amu/cc]'),/nodata  ELSE plot,radius[disk],alog10(g[disk].dens*dens_convert),yrange = [-1.99,3],xrange = xrange,psym = 3,ystyle = 1,/nodata
    oplot,radius[disk],alog10(g[disk].dens*dens_convert*frac),psym = 3,color = obscolor
    loadct,ctables[i]
    oplot,[r25[i],r25[i]],[-4,4],thick = 4,color = colors[i];,linestyle = linestyle[i],color = obscolor
    IF i eq mid THEN ytit_string = textoidl('log \Sigma_{gas} [amu/cm^2]') ELSE ytit_string = ""
    IF i eq n-1 THEN AXIS,YAXIS = 1,YRANGE = [19,24],ytitle = ytit_string,/save,ystyle = 1 ELSE  AXIS,YAXIS = 1,YRANGE = [19.001,24],ytitle = ytit_string,/save,ystyle = 1
    oplot,radiibins/cm_per_kpc,alog10(sd[*,i]),psym = -1*symcat(symbols[i]),thick = 4,color = colors[i];,color = obscolor
;    IF i eq n-1 THEN AXIS,YAXIS = 1,YRANGE = [-1,4],ytitle = ytit_string,/save,ystyle = 1 ELSE  AXIS,YAXIS = 1,YRANGE = [-1.001,4],ytitle = ytit_string,/save,ystyle = 1
;        oplot,radiibins/cm_per_kpc,alog10(sd[*,i]*convert2msol),psym = -1*symcat(symbols[i]),thick = 4;colors[i],color = obscolor
    multiplot,/doyaxis   
ENDFOR
    multiplot,/doyaxis 
multiplot,/reset
formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
IF KEYWORD_SET(color) then loadct,39
END

PRO densRadius_Master,outplot = outplot
prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
tfiles = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
dKpcUnit = [25000.,25000.]
dMsolUnit = [ 2.310e15, 2.310e15]
formatplot,outplot = outplot

densRadius,prefix + tfiles,dKpcUnit,dMsolunit,maxdistance = 8,outplot = outplot

END
