; this is following haynes et al. 1999 section 3.1: a robust width algorithm.
;
; w_21: full width across profile at 50% of each peak:
;       - fit a polynomial between 15% and 85% of respective horn
;       - default straight line. optional 2nd degree.

;incl should be in degrees
;pa should be in degrees
function get_width_fit, vaxis, spectrum, width, v1, v2, widthval1, widthval2,$
                        doplot=doplot, gal=gal, incl=incl, pa=pa, $
                        interactive=interactive,liveplot = liveplot


; get the values where the spectrum is 15% and 85% of its respective
; peak.
;
; a schematic representation of the profile with my variables
; could be:
;
;   |    peak1 --> /\                  /\ <-- peak2
;   |             /  \                /  \
; f |  v1_85 --> /    \              /    \ <-- v2_85
; l |           /      \            /      \ 
; u |          /        \          /        \
; x |         /    o     \________/    o     \
;   | v1 --> /________________________________\ <-- v2
;   |       /              w_50                \
;   |      /                                    \
;   |     /<-- v1_15                    v2_15 -->\
;   |    /                                        \ 
;   ----------------------------------------------------
;                         velocity
;
; 
; We will be fitting a straight line between v1_15 and v1_85 (and a
; second between (v2_85 and v2_15), and using those lines to find
; the location of 50% of the peak on either side.
;

w_15 = get_width_peak(vaxis, spectrum, 0.15, v1=v1_15, v2=v2_15, $
                       peak1=peak1, peak2=peak2, $
                       widthval1=widthval1_15, widthval2=widthval2_15)

w_85 = get_width_peak(vaxis, spectrum, 0.85, v1=v1_85, v2=v2_85, $
                       peak1=peak1, peak2=peak2, $
                       widthval1=widthval1_85, widthval2=widthval2_85)




done = 0
while done le 0 do begin

;**************************** FIRST peak *******************************

    ; this happens if we haven't gone through the loop once already.
    if done eq 0 then begin
        if v1_85 gt v1_15 then ind_1 = where(vaxis ge v1_15 and vaxis le v1_85) $
        else ind_1 = where(vaxis lt v1_85 and spectrum gt peak1 * 0.15)

        ; make sure it found something - if not pick the closest value
        if ind_1[0] eq -1 then begin
            blah = min(abs(vaxis - mean([v1_85, v1_15])), ind_tmp)
            ind_1 = ind_tmp
        endif
    endif
    
    ; this only happens if we've gone through once, determined
    ; the profile to be bad, and then wanted to fit interactively.
    ; this will involve clicking and reading the cursor position.
    if ( (done eq -1) or (done eq -3) ) then begin

        print, "Redo first [left] side of the profile."
        print, "Click on bottom and top region to use for the fits."
        print, "Right mouse click to keep the fit."

        cursor, x1, y1, /down, /data

        ; check for a right click
        if (!mouse.button ne 4) then begin
            ; get the second side of the box
            cursor, x2, y2, /down, /data

            ; overplot the lines
            oplot, [x1, x1], [-1, 1e15], linestyle=2
            oplot, [x2, x2], [-1, 1e15], linestyle=2

            ; make sure the user clicked in the right order - fix if not.
            if x1 gt x2 then begin
                tmp = x1
                x1 = x2
                x2 = tmp
            endif
            ; get the indices in between the cursor positions
            ind_1 = where(vaxis gt x1 and vaxis lt x2)
            
            done = 1
        endif else begin
            if done eq -3 then done = -2
            if done eq -1 then done = 1
        endelse

    endif
    

    ; add in adjacent points if there are only two points to fit to
    if n_elements(ind_1) eq 1 then begin 
        tmp = fltarr(3)
        tmp[0] = ind_1[0] - 1
        tmp[1] = ind_1[0]
        tmp[2] = ind_1[0] + 1
        
        ind_1 = tmp
    endif
    if n_elements(ind_1) eq 2 then begin
        tmp = fltarr(4)
        tmp[0] = ind_1[0] - 1
        tmp[1] = ind_1[0]
        tmp[2] = ind_1[1]
        tmp[3] = ind_1[1] + 1
        
        ind_1 = tmp
    endif


    ;***************************** SECOND peak *****************************
    
    ; figure out where velocity is between 85% and 15% of peak
    if (done eq 0) then begin
        if v2_85 lt v2_15 then ind_2 = where(vaxis le v2_15 and vaxis gt v2_85) $
        else ind_2 = where(vaxis gt v2_85 and spectrum gt peak2 * 0.15)

        if ind_2[0] eq -1 then begin
            blah = min(abs(vaxis - mean([v2_85, v2_15])), ind_tmp)
            ind_2 = ind_tmp
        endif
 
   endif

    ; as before this only happens with a bad fit in interactive mode.
    if ( (done eq -2) or (done eq -3) ) then begin

        print, "Redo second [right] side of the profile."
        print, "Click on bottom and top region to use for the fits."
        print, "Right mouse click to keep the fit."

        cursor, x1, y1, /down, /data

        ; check for a right click
        if (!mouse.button ne 4) then begin
            ; get the second side of the box
            cursor, x2, y2, /down, /data

            ; overplot the lines
            oplot, [x1, x1], [-1, 1e15], linestyle=2
            oplot, [x2, x2], [-1, 1e15], linestyle=2

            ; make sure the user clicked in the right order - fix if not.
            if x1 gt x2 then begin
                tmp = x1
                x1 = x2
                x2 = tmp
            endif
            ; get the indices in between the cursor positions
            ind_2 = where(vaxis gt x1 and vaxis lt x2)
            done = 1
        endif else done = 1

    endif   ; end of interactive mode.

    
    ; add in adjacent points if there are only one or two points to fit
    if n_elements(ind_2) eq 1 then begin 
        tmp = fltarr(3)
        tmp[0] = ind_2[0] - 1
        tmp[1] = ind_2[0]
        tmp[2] = ind_2[0] + 1
        
        ind_2 = tmp
    endif
    if n_elements(ind_2) eq 2 then begin
        tmp = fltarr(4)
        tmp[0] = ind_2[0] - 1
        tmp[1] = ind_2[0]
        tmp[2] = ind_2[1]
        tmp[3] = ind_2[1] + 1
        
        ind_2 = tmp
        
    endif   
    

    ; ********************** FIT the polynomials to both peaks *********************

    ; now do the fits for both.    
    x_1 = vaxis[ind_1]
    y_1 = spectrum[ind_1]
    
    fit_1 = poly_fit(x_1, y_1, 1, yfit=yf_1, chisq=chisq1)

    x_2 = vaxis[ind_2]
    y_2 = spectrum[ind_2]
    fit_2 = poly_fit(x_2, y_2, 1, yfit=yf_2, chisq=chisq2)



    ; now re-arrange on each side to find the width.
    widthval1 = peak1 * 0.5
    widthval2 = peak2 *0.5
    v1 = (peak1 * 0.5 - fit_1[0]) / fit_1[1]
    v2 = (peak2 * 0.5 - fit_2[0]) / fit_2[1]

    
    ; plot to make sure it worked.
    if keyword_set(liveplot) then begin
;        !p.multi=[0,1,2]
        set_plot,'x'
        
        plot,vaxis,spectrum, $
          xtitle='l.o.s. velocity [km / s]', $
          ytitle = 'mass [solar masses]', $
          title = gal + ': i=' + strtrim(round(incl),1) + $
          ', pa=' + strtrim(round(pa),1),xrange = [-200,200]
        
        oplot, x_1, yf_1, color=100
        oplot, x_2, yf_2, color=100
        oplot, [v1, v2], [peak1, peak2] * 0.5, color=250
        oplot, [v1, v2], [peak1, peak2] * 0.5, psym=1, color=250
        oplot, [v1_15, v1_85, v2_15, v2_85], $
          [widthval1_15, widthval1_85, widthval2_15, widthval2_85], $
          psym=1, color=100
; don't actually plot the derivative        
;        cumul = total(spectrum, /cumulative)
;        cumul = cumul / max(cumul)
;        plot, vaxis, cumul, xtitle='velocity [km/s]', ytitle='cumul spectrum'
;        oplot, [v1, v1], [-1, 1e15], color=250
;        oplot, [v2, v2], [-1, 1e15], color=250
        
;        oplot,vaxis[ind_1], cumul[ind_1], color=100
;        oplot, vaxis[ind_2], cumul[ind_2], color=100
        
        
    endif
    
;    print, 'Width is: ' + string(abs(v1-v2)) + ' km/s'
    
    
    ; check the fits if interactive is set.
    ; set done eq 0 to make sure we only go through this loop once!!
    test1 = abs(total((spectrum[ind_1]-yf_1)/spectrum[ind_1])) 
    test2 = abs(total((spectrum[ind_2]-yf_2)/spectrum[ind_2])) 
;            print, '1:',test1
;            print, '2:',test2

    if keyword_set(interactive) and done eq 0 then begin
        if (test1 gt 1.5) and (test2 gt 1.5) then begin
            print, '1: fail'
            print, '2: fail'
            done = -3

            
        endif
        if (test1 gt 1.5) then begin 
            print, '1: fail'
            done = -1

            
        endif
        if (test2 gt 1.5) then begin
            print, '2: fail'
            done = -2

            
        endif
    endif
    if done eq 0 then done = 1
    


endwhile

; make a hard copy of the plots if requested
if keyword_set(doplot) then begin
    !p.multi=0
    set_plot,'ps'
    loadct,39,/silent
    device,filename = './' + gal + '.' + strtrim(round(incl),1) + $
      '.' + strtrim(round(pa*180/!pi),1) + '.ps', /color
    plot,vaxis,spectrum, $
      xtitle='l.o.s. velocity [km / s]', $
      ytitle = 'mass [solar masses]', $
      title = gal + ': i=' + strtrim(round(incl),1) + $
      ', pa=' + strtrim(round(pa),1),xrange = [-200,200]
    
    oplot, x_1, yf_1, color=60
    oplot, x_2, yf_2, color=60
    oplot, [v1, v2], [peak1, peak2] * 0.5, color=250
    oplot, [v1, v2], [peak1, peak2] * 0.5, psym=1, color=250
    device,/close
endif

return, abs(v1-v2)
    
end
