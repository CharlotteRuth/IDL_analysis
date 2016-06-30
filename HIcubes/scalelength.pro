FUNCTION scalelength, star, outplot = outplot,verbose = verbose,rmin = rmin,rmax = rmax

IF keyword_set(verbose) THEN BEGIN
    loadct,0
    !p.multi = 0
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
        !x.charsize=1.5         ;2.25
        !y.charsize=1.5         ;2.25
        !X.MARGIN = [12,3]
        !Y.MARGIN = [6,2]
    ENDIF ELSE BEGIN
        !P.CHARTHICK=1.5
        !X.THICK=1.5
        !Y.THICK=1.5
        !p.charsize=1.0
        !x.charsize=1.5
        !y.charsize=1.5  
        !X.MARGIN = [12,3]
        !Y.MARGIN = [6,2]
        set_plot,'x'
        nbins=100.0
        linestyles = [0,2]
    ENDELSE
ENDIF

IF NOT KEYWORD_SET(rmin) THEN rmin = 0
IF NOT KEYWORD_SET(rmax) THEN rmax = 15

range = rmin - rmax
prof = prof(star,'star',MAX(star.tform),nbins = nbins, rmax = rmax)
;prof_stellar[0,iFile,*] = prof.rbins
;prof_stellar[1,iFile,*] = prof.rho 

massform = MAX(star.mass)
fitpar =  [0,0,0,0,0]
ind = where(prof.rbins gt rmin and prof.rbins lt rmax)   
prof.rho = prof.rho/1e6 ;Now in pc



;dblexpfit, prof.rbins[ind], prof.rho[ind], $
;  prof.rho[ind]/sqrt(prof.num_in_bin[ind]), fitpar, red_chisq=exp2chi, /no_guess, do_plot = 1
fitpar = poly_fit(prof.rbins[ind], alog(prof.rho[ind]),1)

IF KEYWORD_SET(verbose) THEN BEGIN
    plot,prof.rbins,prof.rho,/ylog,xrange = [0.1,rmax]
    oplot,prof.rbins,exp(fitpar[0])*exp(prof.rbins*fitpar[1]),linestyle = 2
;oplot,prof.rbins,fitpar[3]*exp(-1.0*prof.rbins/fitpar[4]),linestyle = 2
ENDIF
return,-1.0/fitpar[1]
END
