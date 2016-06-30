PRO LGanalysis
stepst = ['504','512']
files = '/astro/net/scratch2/christensen/LG/gas/steps/LG.'+stepst+'.dir/LG.'+stepst
outplot = '/astro/net/nbody1/christensen/Schmidtlaw/data/LG.504.512.guo.ps'
psym = [4,5]
color = [50,240]
mulitHaloGuo, files, outplot = outplot, color = color, psym = psym 
END

PRO h516analysis
dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/'
;step = '00168'
;files = ['h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.'+step+'.dir/h516.cosmo25cmb.1536g1MBWK.'+step, $
;         'h516.cosmo25cmb.1536g2HBWK_nocool/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step, $
;         'h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step]
;step = '00512'
;files = ['h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.'+step+'.dir/h516.cosmo25cmb.1536g1MBWK.'+step, $
;         'h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step]
step = '00168'
files = ['h516.cosmo25cmb.1536g2HBWK_nocool/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step, $
         'h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step]
key = ['H2, w/o cooling', 'H2 w/ cooling']
files = dir + files
multiHaloGuo,files,color = [50,240],psym = [4,4],key = key

END

PRO mulitHaloGuo, files, outplot = outplot, color = color, psym = psym, key = key 
loadct,39
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23

IF KEYWORD_SET(outplot) THEN BEGIN
   set_plot, 'ps'
   !P.THICK=1.5                 ;4
   !P.CHARTHICK=1.5             ;4
   !X.THICK=1.5                 ;4
   !Y.THICK=1.5                 ;4
   !p.charsize=1.2              ;1.8
   !p.font=0
;   device, filename=outbase+'LG.'+step[0]+'.guo.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
   device, filename=outplot,/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
ENDIF ELSE BEGIN
    set_plot,'x'
    window,1
ENDELSE

IF NOT KEYWORD_SET(psym) THEN psym = 2
IF N_ELEMENTS(psym) eq 1 THEN psym = fltarr(N_ELEMENTS(step)) + psym 
IF NOT KEYWORD_SET(color) THEN color = 240
IF N_ELEMENTS(color) eq 1 THEN color = fltarr(N_ELEMENTS(step)) + color 

x = findgen(700)/100 + 9
y = alog10(0.1257*((10^x/10^(11.36))^(-0.9147) + (10^x/10^(11.36))^(0.2485))^(-2.574)*10^x) ;guo 09, (3)
plot,x,y,xstyle = 1, ystyle = 1,xrange = [9.5,13],yrange = [6,11],xtitle = textoidl('log(M_{halo}[M')+sunsymbol()+'])',ytitle = textoidl('log(M_{*}[M')+sunsymbol()+'])'

FOR j = 0, N_ELEMENTS(files) -1 DO BEGIN
    file = files[j]
    nhalos = 100

    guo_mass = fltarr(2,nhalos)
    FMT = 'X,X,X,X,X,F,X,F,F,F,X,X,X,X,X,X,X,X,X,X,A'
    readcol,file + '.amiga.stat',F=FMT,vir,gas,star,dark,contam,/silent 
    FOR i = 0, nhalos - 1 DO BEGIN
        guo_mass[*,i] = [dark[i],star[i]]
    ENDFOR
    IF NOT KEYWORD_SET(psym) THEN psym = 2
    IF N_ELEMENTS(psym) eq 1 THEN psymi = psym ELSE psymi = psym[j] 
    IF NOT KEYWORD_SET(psym) THEN psym = 2
    IF N_ELEMENTS(psym) eq 1 THEN psymi = psym ELSE psymi = psym[j] 

    oplot,alog10(guo_mass[0,*]),alog10(guo_mass[1,*]),psym = psym[j],color = color[j]
ENDFOR
legend,key,psym = psym,color = color,/bottom,/right
IF KEYWORD_SET(outplot) THEN device,/close

;IF KEYWORD_SET(outplot) THEN BEGIN
;   set_plot, 'ps'
;   !P.THICK=1.5                 ;4
;   !P.CHARTHICK=1.5             ;4
;   !X.THICK=1.5                 ;4
;   !Y.THICK=1.5                 ;4
;   !p.charsize=1.2              ;1.8
;   !p.font=0
;;   device, filename=outbase+'LG.'+step[0]+'.guo.ps',/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
;   device, filename=outplot,/COLOR,bits_per_pixel= 8,/times,ysize=3.5,xsize=5,/inch
;ENDIF ELSE BEGIN
;    set_plot,'x'
;    window,1
;ENDELSE

;plot,[0,0],[1,1],xrange=[1e-10,1e5],yrange=[1e1,1e8],/ylog,/xlog
;massunit = 9.9669000e+16
;lengthunit = 87671.200
;rhounit=massunit * gm_per_msol * H_per_gm/lengthunit^3/cm_per_kpc^3
;FOR j = 0, N_ELEMENTS(step) -1 DO BEGIN
;    stepst = step[j]
;    dir = '/astro/net/scratch2/christensen/LG/'
;    file = stepst
;    nhalos = 30

;    rtipsy,dir+'/'+file,h,g,d,s
;    IF NOT KEYWORD_SET(psym) THEN psym = 2
;    IF N_ELEMENTS(psym) eq 1 THEN psymi = psym ELSE psymi = psym[j] 
;    IF NOT KEYWORD_SET(psym) THEN psym = 2
;    IF N_ELEMENTS(psym) eq 1 THEN psymi = psym ELSE psymi = psym[j] 

;    oplot,g.dens*rhounit,g.tempg,psym = 3,color = color[j]
;ENDFOR
;legend,'LG.'+step,psym = psym,color = color,/bottom,/right
;IF KEYWORD_SET(outplot) THEN device,/close ELSE stop

END
