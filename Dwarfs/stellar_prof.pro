pro stellar_prof,infile,outfile,densmin,densmax,SECONDFILE = secondfile,TITLE = TITLE
;This plots a series of stellar profiles showing the evolution in time
;stellar_prof,'../1E5R/10M/o10M_1.00300','10M_sprof',1e2,1e8,secondfile=
;'../1E5R/10M_original/o10M_1.00300',title = '1e10 Solar Mass' 

msol = 2.362e5
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
timeunit = 1e9
nframes = 20
maxtime = 3e9
nbins = 25
rmax = 10
rmin = 0

rtipsy,infile,h,g,d,s
IF (KEYWORD_SET(secondfile)) THEN rtipsy,secondfile,h2,g2,d2,s2
IF NOT (KEYWORD_SET(TITLE)) THEN title = infile

time = (findgen(nframes)+1.)*maxtime/nframes
files = sindgen(nframes)

FOR timect = 0,nframes-1 DO BEGIN
    ind = WHERE(s.tform GT 0.0 AND s.tform*timeunit LT time[timect])
    s_good = s[ind]
    s_prof = prof(s_good, 'star', h.time, nbins = nbins, rmax = rmax)
;    fitpar =  [0,0,0,0,0]    
    s_prof.rho = s_prof.rho * msol ; Msol/kpc^2

    IF (KEYWORD_SET(secondfile)) THEN BEGIN
        ind2 = WHERE(s2.tform GT 0.0 AND s2.tform*timeunit LT time[timect])
        s2_good = s2[ind2]
        s2_prof = prof(s2_good, 'star', h2.time, nbins = nbins, rmax = rmax)
;    fitpar =  [0,0,0,0,0]    
        s2_prof.rho = s2_prof.rho * msol ; Msol/kpc^2
    ENDIF

    plot,s_prof.rbins,s_prof.rho,/ylog,title = title + STRTRIM(time[timect]),xtitle = 'Radius',ytitle = 'Stellar Density',yrange = [densmin,densmax]
    IF (KEYWORD_SET(secondfile)) THEN BEGIN
        oplot,s2_prof.rbins,s2_prof.rho,linestyle = 1 
        legend,['Higher SF Threshold','Original'],linestyle=[0,1]
    ENDIF   

    set_plot,'ps'
    device,filename=outfile+STRTRIM(files[timect])+'.eps',/color,bits_per_pixel=8 
    plot,s_prof.rbins,s_prof.rho,/ylog,title = title + STRTRIM(time[timect]),xtitle = 'Radius',ytitle = 'Stellar Density',yrange = [densmin,densmax]
    IF (KEYWORD_SET(secondfile)) THEN BEGIN
        oplot,s2_prof.rbins,s2_prof.rho,linestyle = 1
        legend,['Higher SF Threshold','Original'],linestyle=[0,1]
    ENDIF
    device,/close
    set_plot,'x'
ENDFOR

END
