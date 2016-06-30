;Based on Ferah's lineflux program
;Returns the Halpha flux in ergs/s
;SFR=7.9e-42* L(halpha).
function Halpha, filename,outplot = outplot,verbose = verbose,extno_int = extno_int
;filename is the name of the mcrx file
IF NOT KEYWORD_SET(extno_int) THEN extno_int = 32 ELSE extno_int = extno_int
cube=mrdfits(filename, extno_int) ;changes depending on mcrx file; Should be integrated quantities
x=cube.L_LAMBDA
lambda=cube.LAMBDA
w=where(cube.LAMBDA gt 6.e-7 and lambda lt 7.e-7);trim the spectrum
x=x[w]
lambda=lambda[w]
rl = where(lambda gt 6.54e-7 and lambda lt 6.58e-7,nl)
rc = where((lambda gt 6.e-7 and lambda le 6.5e-7) or $ ;define continuum
             (lambda gt 6.65e-7 and lambda lt 7.e-7))
;con = mean(x[rc])
con = SPLINE(lambda[rc],x[rc]*1d-40,6.56e-7)
con = con*1d40
fl = $
    tsum(lambda,x-con,min(rl),max(rl)) ;measure lineflux
print, 'continuum and flux equal',con, fl*1.e7; convert to ergs/s

IF (KEYWORD_SET(verbose)) THEN BEGIN
;    loadct,39
    !p.multi = 0
    !Y.STYLE = 1
    !X.STYLE = 1
    !P.THICK = 3.5
    IF KEYWORD_SET(outplot) THEN BEGIN
        !P.CHARTHICK=4
        !X.THICK=4
        !Y.THICK=4
        !p.charsize=1.0
        !x.charsize=1.5         ;2.25
        !y.charsize=1.5         ;2.25
        !X.MARGIN = [12,3]
        !Y.MARGIN = [6,2]
        black  = 0
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
    ENDELSE
    if (KEYWORD_SET(outplot)) then begin
        set_plot,'ps'
        device,filename = outplot+'halpha.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2
    endif else begin
        set_plot,'x'
        window,1
    endelse
    plot, cube.LAMBDA, cube.L_LAMBDA, xrange=[6.e-7, 7.e-7], thick=4, xtit='Wavelength(m)', ytit='Flux (W/m)'
    oplot, lambda, x,color = 240
    oplot, lambda, lambda*0 + con, linestyle = 1
    if (KEYWORD_SET(outplot)) then device,/close
    stop
ENDIF    
return,fl*1.e7
END
