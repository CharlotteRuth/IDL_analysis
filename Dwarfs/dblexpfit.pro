; attempts to fit two exponentials to a radial density profile


;  *** USAGE ***
; 
;  The procedure dblexpfit takes 7 parameters:
;
;  x, y
;  err - error for each point, must be an array of same size as x,y

; array p of initial guesses, with the components:

; p[0] = break radius (not required for the fit, but returned for practicality)
; p[1] = sfc brightness of inner expo
; p[2] = scale radius of inner expo
; p[3] = sfc brightness of outer expo
; p[4] = scale radius of outer expo
;
; if no initial guess is given, a crude linear fit is done and tweaked
; slightly to get initial guesses - this works if the second
;                                   exponential is steeper than the first


FUNCTION dblexp, p, X=x, Y=y, ERR=err
br = (alog(p[1]) - alog(p[3]))/(1./p[2]-1./p[4])

res = fltarr(n_elements(x))

for i=0,n_elements(x)-1 do begin
    if x[i] lt br then res[i] = p[1]*exp(-x[i]/p[2]) $
      else res[i] = p[3]*exp(-x[i]/p[4])
endfor

return, (y - res)/err
;return, y - res
end



pro dblexpfit, x, y, err, p0, DO_PLOT = do_plot, CHISQ = chisq, RED_CHISQ=red_chisq, NO_GUESS = no_guess

; if no guesses are provided, do a really stupid initial guess
if(n_params() lt 4 or keyword_set(NO_GUESS)) then begin
    nelem = n_elements(y)

    logx = alog(x)
    logy = alog(y)
    slope = (logy[nelem-1]-logy[0])/(x[nelem-1]-x[0])
    b = logy[nelem/2]-slope*x[nelem/2]
    ; the rest of the stuff assumes slope > 0
    if(slope lt 0.) then slope = -slope
    print, slope, 1./slope, exp(b)
    r_br = total(x)/n_elements(x)    
    r1 = 1./slope
    sig1 = exp(b)
    r2 = 1./slope/2.
    sig2 = exp(alog(sig1)-r_br*(1./r1 - 1./r2))
    p0=[r_br, sig1, r1, sig2, r2]
    ;p0=[5.,exp(b), 1./slope, exp(b)*100., 1./slope/2.]
endif
    
fa = {X:x, Y:y, ERR: err}

; set up par info to keep break radius fixed
parinfo = replicate({fixed:0},5)
parinfo[0].fixed = 1
p = mpfit('dblexp', p0, functargs=fa, $
          BESTNORM = chisq, PERROR = perror, PARINFO = parinfo)
chisq = chisq^2 ; the bestnorm returned is just the residuals, not chisq
; reduced chisq = chisq / (N_points - N_parameters - 1)
red_chisq = chisq/(n_elements(x) - n_elements(p) - 1.)

p[0] = (alog(p[1]) - alog(p[3]))/(1./p[2]-1./p[4])
print, 'r_br  sig_in   r_in   sig_out   r_out'
print, p 
print, 'r_br/r_in = ', p[0]/p[2]
print, 'chisq: ', chisq, ' reduced chisq: ', red_chisq

if(keyword_set(do_plot)) then begin

    !p.multi = [0,1,2]

    plot, x, y, psym = 2, /ylog, xrange = [min(x),max(x)], title = 'initial guesses',$
      yrange = [1,1000], ystyle  = 1
    oplot, x, p0[1]*exp(-x/p0[2]), line = 3
    oplot, x, p0[3]*exp(-x/p0[4]), line = 2
    oploterr, x, y, err
    legend, /right, line=[3,2], ['inner', 'outer']

    plot, x, y, psym = 2, xrange = [min(x),max(x)], /ylog, $
      title = 'fitted reduced chisq = ' + string(red_chisq), yrange = [1,1000], ystyle = 1
    oplot, x, p[1]*exp(-x/p[2]), line = 3
    oplot, x, p[3]*exp(-x/p[4]), line = 2
    oploterr, x, y, err    
    legend, /right, line=[3,2], ['inner', 'outer']
    
    

    !p.multi = 0

endif

p0 = p

return

end
