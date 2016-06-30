PRO color_scaleLength,simdata,outplot = outplot,color = color,symbols = symbols, thicks = thicks,key = key
n = N_ELEMENTS(simdata)
!p.multi = 0
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    black  = 0
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    black = 255
    fgcolor = 255
    bgcolor = 0
ENDELSE
if (KEYWORD_SET(outplot)) then begin
    set_plot,'ps'
    device,filename = outplot+'Magscale.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
endif else begin
    set_plot,'x'
    window,0,xsize = 712,ysize = 392
endelse
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  color  = (findgen(N_ELEMENTS(broadband)) + 1)*240/N_ELEMENTS(broadband) else color = color
    IF NOT KEYWORD_SET(thicks) THEN thicks = fltarr(n) + 2
    IF NOT KEYWORD_SET(symbols) THEN symbols = fltarr(n) + 4
    obscolor = fgcolor
    obssym = 4
ENDIF ELSE BEGIN
    loadct,0    
    color = fltarr(n) + fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT KEYWORD_SET(thicks) THEN thicks = (fltarr(n) + 1)
    IF NOT KEYWORD_SET(symbols) THEN symbols = (findgen(n)+2)*2
    obscolor = 150
    obssym = 2
ENDELSE

datafile = '~/code/Sunrise_analysis/vanZee_short.txt'
readcol,datafile,name1,name2,D25,d_25,I,PA,m_B,m_B_error,A_B,B_V, B_V_error,U_B,U_B_error,D,MB,Sigma_25,mu0B,mu0cB,alpha_kpc,alpha_arcsec,Halpha,Halpha_error,LF,format='(A,A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)' 

!p.multi = [0];,1,2]
plot,MB,alog10(alpha_kpc),psym = obssym,xtitle = textoidl('M_B'),ytitle = textoidl('log \alpha (kpc)')
oplot,MB,alog10(alpha_kpc),psym = obssym,color = obscolor
FOR i=0,n -1 DO BEGIN
;    oplot,[simdata[i].B,simdata[i].B],alog10([simdata[i].alpha,simdata[i].alpha]),psym = symbols[i],color = color[i]
    oplot,[simdata[i].B,simdata[i].B],alog10([simdata[i].alpha_gal,simdata[i].alpha_gal]),psym = symbols[i],color = color[i]
ENDFOR
IF KEYWORD_SET(key) THEN legend,['van Zee `00',key],color = [obscolor,color],psym = [obssym,symbols],/right


IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot+'Magcolor.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
ENDIF ELSE     window,1,xsize = 712,ysize = 392

plot,MB,B_V,psym = obssym,xtitle = textoidl('M_B'),ytitle = textoidl('B-V'),yrange=[0.2,0.8]
oplot,MB,B_V,psym = obssym,color = obscolor
FOR i=0,n -1 DO BEGIN
    oplot,[simdata[i].B,simdata[i].B],[simdata[i].B - simdata[i].V, simdata[i].B - simdata[i].V],psym = symbols[i],color = color[i]
ENDFOR
IF KEYWORD_SET(key) THEN legend,['van Zee `00',key],color = [obscolor,color],psym = [obssym,symbols],/right

IF (KEYWORD_SET(outplot)) THEN BEGIN
    device,/close
    device,filename = outplot+'Scalecolor.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2
ENDIF ELSE     window,3,xsize = 712,ysize = 392

plot,B_V,alog10(alpha_kpc),psym = obssym,xtitle = textoidl('M_B'),ytitle = textoidl('B-V');,yrange=[0.2,0.8]
oplot,B_V,alog10(alpha_kpc),psym = obssym,color = obscolor
FOR i=0,n -1 DO BEGIN
    oplot,[simdata[i].B - simdata[i].V, simdata[i].B - simdata[i].V],alog10([simdata[i].alpha_gal,simdata[i].alpha_gal]),psym = symbols[i],color = color[i]
ENDFOR
IF KEYWORD_SET(key) THEN legend,['van Zee `00',key],color = [obscolor,color],psym = [obssym,symbols],/right

if (KEYWORD_SET(outplot)) then device,/close else stop
!p.multi = 0
END



PRO color_scaleLength_master
;gh516_H = {B:-16.5029574262822, V:-16.8630676638494,alpha:1.1024871,alpha_gal:1.48329,alphaV:1.1310378,alpha_galV:1.51750} ;2304
gh516_H = {B:-15.6996847198888, V:-16.0784155581285, alpha:0.81438442,alpha_gal:1.02964, alphaV:0.00000001,alpha_galV:0.000001} ;h516.cosmo25cmb.3072g14HBWK
gh516_M = {B:-15.2652157189731, V:-15.7550592733375, alpha:0.79430061,alpha_gal:0.770795,alphaV:0.79689074,alpha_galV:0.777515}
simdata = [gh516_M,gh516_H]
key = ['DnoH2','DH2']
outplot = '~/plots/h516.cosmo25cmb.paper_'
color_scaleLength,simdata,key = key,outplot = outplot
END
