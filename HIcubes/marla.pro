pro marla, color=color

; read in and select the data.
; selection critera are marla's.
datadir = '/astro/users/christensen/code/HIcubes/'
sfi = mrdfits(datadir + 'SFI_data.fits',1)
q=where(sfi.W20_err/sfi.W20 le 0.01 and sfi.absmag_err le 0.2 and $
        sfi.sini ge 0.7,n)
sfi = sfi[q]
verh = mrdfits(datadir + 'verh_data.fits',1)
q=where(verh.w20_err lt 10)
verh=verh[q]
matt = mrdfits(datadir + 'matthews_data.fits',1)
q=where(matt.sini ge 0.6)
matt=matt[q]
mcg = mrdfits(datadir + 'mcg_data.fits',1)

; read in simulation data.
readcol, datadir + 'marla_data.txt', galnum, mw20, mw20err, mi, mierr,/silent
readcol, datadir + 'tf.txt', w20, i,/silent


; create structures to store the data.
simdata = replicate( $
                     create_struct('w20', 0.0, $
                                   'w20_err', 0.0, $
                                   'i', 0.0, $
                                   'i_err', 0.0), $
                     n_elements(w20))


h258 = replicate( $
                     create_struct('w20', 0.0, $
                                   'w20_err', 0.0,$
                                   'i', 0.0, $
                                   'i_err', 0.0), $
                     4)

marladata = replicate( $
                       create_struct('w20', 0.0, $
                                     'w20_err', 0.0, $
                                     'i', 0.0, $
                                     'i_err', 0.0), $
                       n_elements(mw20))

sfidata = replicate( $
                     create_struct('w20', 0.0, $
                                   'w20_err', 0.0,$
                                   'i', 0.0, $
                                   'i_err', 0.0), $
                     n_elements(sfi))

verhdata = replicate($
                      create_struct('w20', 0.0, $
                                    'w20_err', 0.0,$
                                    'i', 0.0, $
                                    'i_err', 0.0), $
                      n_elements(verh))

mattdata = replicate( $
                      create_struct('w20', 0.0, $
                                    'w20_err', 0.0,$
                                    'i', 0.0, $
                                    'i_err', 0.0), $
                    n_elements(matt))

mcgdata = replicate( $
                     create_struct('w20', 0.0, $
                                   'w20_err', 0.0,$
                                   'i', 0.0, $
                                   'i_err', 0.0), $
                     n_elements(mcg))


; put data in structures. subtract 8 km/s for the turbulence correction.
simdata.w20 = w20 - 8
simdata.i = i
; 0.8, 0.4, 0.2, 0
h258.w20 = [244.6, 246.6, 242.7, 245] - 8
;h258.i = [-22.91, -22.66, -22.58, -22.6]
;h258.i = [-23.10, -22.70, -22.64, -22.49]
h258.i = [-23.0, -22.6, -22.5, -22.4]

marladata.w20 = mw20 - 6
marladata.w20_err = mw20err
marladata.i = mi
marladata.i_err = mierr


sfidata.w20 = sfi.w20/2.
sfidata.w20_err = sfi.w20_err
sfidata.i = sfi.absmagi

verhdata.w20 = verh.w20/2.
verhdata.w20_err = verh.w20_err
verhdata.i = verh.absmag_bvri[3]

mattdata.w20 = matt.w20/2.
mattdata.w20_err = matt.w20_err
mattdata.i = matt.absmagi


mcgdata.w20 = mcg.w20
mcgdata.w20_err = mcg.w20_err
mcgdata.i = mcg.absmagi


; set up color scheme
col = indgen(8)

cmcg = 0
cmatt = 1
csfi = 2
cverh = 3
cmarla = 4
csim = 5
ch258 = 6
ch2582 = 7
if keyword_set(color) then begin
    loadct,39
    col = [50, 80, 120, 160, 180, 250, 150, 0]
endif else begin
    col = [100, 100, 100, 100, 100, 0, 0, 250]
endelse

a = findgen(17) * (!pi * 2/16.)
usersym, 0.2*cos(a), 0.2*sin(a), /fill

; plot the data.

plot, mcgdata.w20, mcgdata.i, /xlog, yrange=[-13.5,-24.1], $
  xrange=[9, 400], /nodata, ystyle=1, xstyle=1, $
  xtitle='W' + textoidl('_{20}') + ' / 2 [km s' + textoidl('^{-1}') + ']', $
  ytitle='M' + textoidl('_I') + ' - 5' + $
  textoidl('log h_{70}'), thick=4
oplot, mcgdata.w20, mcgdata.i, psym=8, color=col[cmcg]
oplot, mattdata.w20, mattdata.i, psym=8, color=col[cmatt]
oplot, sfidata.w20, sfidata.i, psym=8, color=col[csfi]
oplot, verhdata.w20, verhdata.i, psym=8, color=col[cverh]
oplot, marladata.w20, marladata.i, psym=8, color=col[cmarla]
oploterror, marladata.w20, marladata.i, marladata.w20_err, marladata.i_err,psym=3, color=100


; Fits from numbers in paper FAIL.
;w20fitp = [10,25,300,500]
;a = 0.37
;b = -0.094
;ap = -1. * a / b
;bp = 1 / b
;ap = 4.11
;bp = -10
;ifitp = ap + bp* alog10(w20fitp/ 50.) - 20 + 5*alog10(73/70.)
;oplot,w20fitp, ifitp, color=255
;ifit = [-12, -14, -24, -25]
;w20fit = 10.^(a + b * (ifit + 20)) * 50.
;oplot,w20fit,ifit, color=255


; fit from looking at the figure
x = [30,300]
y = [-14.5, -24]
coeff = line(alog10(x), y)
x = [10, 500]
y = coeff[0] + coeff[1] * alog10(x)
oplot,x,y, thick=4


; make a star.
usersym, 1.5*[-1, -0.33, 0, 0.33, 1, 0.33, 0, -0.33, -1], $
         1.5*[0, 0.33, 1, 0.33, 0, -0.33, -1, -0.33, 0], /fill

; plot sim data with triangles.
symbols,22,0.4
oplot, simdata.w20, simdata.i, psym=8, color=col[csim]

; plot h258 in a different symbol
a = findgen(17) * (!pi * 2/16.)
usersym, 0.5*cos(a), 0.5*sin(a), /fill
symbols, 2, 0.8
oplot, [(h258.w20)[3], 0.0001], [(h258.i)[3], 0.0001], psym=8, color=col[csim]


legend, ['  Geha et al. (2006)','  Simulated sample','  Merger remnant'], $
  psym=[8,8,8],  symsize=[0.2, 0.4, 1], $
  position=6, charsize=0.7, $
  color=[col[csim], col[csim], col[ch258]]



             
;  psym=[symcat(16), symcat(46), symcat(16)]

; set up inset and corresponding box
x1 = 220.
x2 = 270.
y1 = -22.3
y2 = -23.2
factor = 4
xmin = alog10(9.)
xmax = alog10(400.)
ymin = -13.5
ymax = -24.1
xr1 = (alog10(x1) - xmin) / (xmax - xmin)
xr2 = (alog10(x2) - xmin) / (xmax - xmin)
yr1 = (y1 - ymin) / (ymax - ymin)
yr2 = (y2 - ymin) / (ymax - ymin)
dx = xr2 - xr1
dy = yr2 - yr1
xb1 = 0.35
yb1 = 0.55
xb2 = xb1 + factor*dx
yb2 = yb1 + factor*dy
pos = [xb1, yb1, xb2, yb2]


; plot box and inset.
oplot, [x1, x1, x2, x2, x1], [y1, y2, y2, y1, y1], thick=3
; these lines are hardcoded
oplot, [57, x2], [-19, y1]
oplot, [56, x2], [-23.55, y2-0.03]
;stop
plot,h258.w20, h258.i, psym=3, position=pos, /noerase, /xlog,$
  xrange=[x1, x2], yrange=[y1, y2], xstyle=1,ystyle=1, charsize=0.5, $
  /nodata,yticks=3

; replot the data in the new box.
a = findgen(17) * (!pi * 2/16.)
usersym, 0.2*cos(a), 0.2*sin(a), /fill
oplot, mcgdata.w20, mcgdata.i, psym=8, color=col[cmcg]
oplot, mattdata.w20, mattdata.i, psym=8, color=col[cmcg]
oplot, sfidata.w20, sfidata.i, psym=8, color=col[cmatt]
oplot, verhdata.w20, verhdata.i, psym=8, color=col[csfi]
oplot, marladata.w20, marladata.i, psym=8, color=col[cmarla]

; plot sim data
a = findgen(17) * (!pi * 2/16.)
usersym, 0.7*cos(a), 0.7*sin(a), /fill
oplot,h258.w20, h258.i, psym=8, color=col[ch258]


oplot,x,y, thick=4

axis, xaxis=0, xtickv=[220,240, 260],xticks=4,/data,charsize=0.5


xyouts, 242, -22.35, 'z = 0', /data, charsize=0.8
xyouts, 239, -22.47, 'z = 0.2', /data, charsize=0.8
xyouts, 241, -22.6, 'z = 0.4', /data, charsize=0.8
xyouts, 235, -22.93, 'z =  0.8', /data, charsize=0.8

; lowest: z=0
;         z=0.2
;         z=0.4
; hi-est: z=0.8



; marla's code:
;all0=create_struct('data', 0, $
;                   'imag', 0., $
;                   'imag_err', 0., $
;                   'sample', -1L)
;all=replicate(all0, n_elements(sfi) + n_elements(verh) + n_elements(matt)+ $
;              n_elements(mcg) + n_elements(data));
;
;;imag= reform(absmag[3,*] - 0.0647 - 0.7177*(absmag[3,*]-absmag[4,*]-0.208), $
;             n_elements(data))
;
;all.vmax = [sfi.W20/2., verh.W20/2.,matt.W20/2,mcg.W20,data.vmax]
;all.vmax_err  = [fltarr(n_elements(sfi))+1.e-3, $
;                 fltarr(n_elements(verh))+1.e-3, $
;                 fltarr(n_elements(matt))+1.e-3, $
;                 fltarr(n_elements(mcg))+1.e-3, $
;                 data.vmax_err]
;all.imag  = [sfi.absmagi, verh.absmag_BVRI[3],matt.absmagi,mcg.absmagi,imag]
;all.imag_err  = [fltarr(n_elements(sfi)), $
;                 fltarr(n_elements(verh)), $
;                fltarr(n_elements(matt)), $
;                fltarr(n_elements(mcg)), $
;                1./sqrt(dodgy.absmag_ivar[3]) ]
;all.mass = [sfi.mass_bary, verh.mass_bary,matt.mass_bary,mcg.mass_bary, $
;            data.mass]
;all.mass_err = [fltarr(n_elements(sfi)), $
;                fltarr(n_elements(verh)), $
;                fltarr(n_elements(matt)), $
;                fltarr(n_elements(mcg)), $
;                data.mass_err]
;all.sample= [replicate(1, n_elements(sfi)), $
;             replicate(2, n_elements(verh)), $
;             replicate(3, n_elements(matt)), $
;             replicate(4, n_elements(mcg)), $
;             replicate(0, n_elements(dodgy))]
;

end
