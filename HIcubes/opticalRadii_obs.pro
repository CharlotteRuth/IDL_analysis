FUNCTION opticalRadii_obs, star, outplot = outplot,verbose = verbose,minr = minr,maxr = maxr

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

IF NOT KEYWORD_SET(minr) THEN minr = 0
IF NOT KEYWORD_SET(maxr) THEN maxr = 20

range = maxr - minr
prof = prof(star,'star',MAX(s.tform),nbins = nbins, rmax = maxr, rmin = minr)
;prof_stellar[0,iFile,*] = prof.rbins
;prof_stellar[1,iFile,*] = prof.rho 

tcurrent = max(star.tform*timeunit)
massform = MAX(star.mass)
fitpar =  [0,0,0,0,0]
ind = where(prof.rbins gt rmin and prof.rbins lt rmax)   
dblexpfit, prof.rbins[ind], prof.rho[ind], $
  prof.rho[ind]/sqrt(prof.num_in_bin[ind]), fitpar, red_chisq=exp2chi, /no_guess

plot,prof.rbins,prof.rho
stop
return,fitpar
END
