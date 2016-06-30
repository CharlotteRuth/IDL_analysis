pro schmidtlaw_res_obs_bigiel,outplot = outplot

res = '0.750000'
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    fgcolor = 0
    bgcolor = 255
ENDIF ELSE BEGIN
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.0
    !y.charsize=1.0  
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
    fgcolor = 255
    bgcolor = 0
ENDELSE
!p.multi = 0

nlevels = 154 
color_hist = findgen(nlevels) + (254 - nlevels)

fileobs = '~/code/HIcubes/bigiel10.dat'
readcol,fileobs,type,name1,name2,logHI,e_logHI,logH2,e_logH2,logSFR,e_logSFR,format = '(A,A,A,F,F,F,F,F,F)'
ind = where(logSFR ne 0)
logSFR = logSFR[ind]
logHI = alog10(10^logHI[ind]/1.36)
logH2 = alog10(10^logH2[ind]/1.36)
logH2[where(logH2 eq 0)] = -6

logH =   alog10((10^(logHI) + 10^(logH2)))
loggas = alog10((10^(logHI) + 10^(logH2))*1.36)
min1 = -0.5
max1 =  2.5
min2 = -4
max2 = -0.5
nbins = 100
bin1 = (max1 - min1)/nbins                      
bin2 = (max2 - min2)/nbins                     
nxbins=FLOOR((max1-min1) / bin1)
nybins=FLOOR((max2-min2) / bin2)
xbins=(FINDGEN(nxbins)*bin1)+min1+bin1/2.0
ybins=(FINDGEN(nybins)*bin2)+min2+bin2/2.0
hist   = HIST_2D(logH, logSFR,min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
histHI = HIST_2D(logHI,logSFR,min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)
histH2 = HIST_2D(logH2,logSFR,min1=min1 + bin1/2,max1=max1 - bin1/2.0,min2=min2,max2=max2 - bin2,bin1 = bin1, bin2 = bin2)

nbins_sim = 50
bin1_sim = (max1 - min1)/nbins_sim                      
bin2_sim = (max2 - min2)/nbins_sim                     
nxbins_sim=FLOOR((max1-min1) / bin1_sim)
nybins_sim=FLOOR((max2-min2) / bin2_sim)
xbins_sim=(FINDGEN(nxbins_sim)*bin1_sim)+min1+bin1_sim/2.0
ybins_sim=(FINDGEN(nybins_sim)*bin2_sim)+min2+bin2_sim/2.0

xsigma = 10.0^(findgen(600)/100 - 1.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4

IF KEYWORD_SET(outplot) THEN device,filename = '~/plots/thesis_resSchmidt_all.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 7,xoffset =  2,yoffset =  2 ELSE window,0,xsize = 1200,ysize = 400
multiplot,[3,1],/square
contour,hist,xbins[0:98],ybins[0:99],xrange = [min1,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill,ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xtitle=textoidl('Log \Sigma')+"!lH!n [M"+sunsymbol()+" pc!u-2!n]", nlevels = nlevels,c_colors = color_hist
multiplot
contour,histHI,xbins[0:98],ybins[0:99],xrange = [-0.49,max1],yrange = [min2,max2], xstyle = 1, ystyle = 1,/fill, xtitle=textoidl('Log \Sigma_{HI} [M')+sunsymbol()+" pc!u-2!n]", nlevels = nlevels,c_colors = color_hist
multiplot
contour,histH2,xbins[0:98],ybins[0:99],xrange = [-0.49,max1],yrange= [min2,max2], xstyle = 1, ystyle = 1,/fill, xtitle=textoidl('Log \Sigma_{H_2} [M')+sunsymbol()+" pc!u-2!n]",nlevels = nlevels,c_colors = color_hist
multiplot,/reset
IF KEYWORD_SET(outplot) THEN device,/close

end
