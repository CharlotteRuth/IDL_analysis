; make some star and gas density plots for run 12M
;example: dens_plots,"masslist.dat",100,1.72043010752688,800


pro dens_plots, list, nbins, sfrthresh, rmax, $
                filename = filename, sfrtitle = sfrtitle, titles = titles, zlimit = zlimit

msol = 2.362e5
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc
readcol, list, format = 'A60,F3,F3', files, rmin, rclip
xrange = [0,14.99]

set_plot,'x'
;set_plot,'ps'
if(keyword_set(filename) eq 0) then filename = 'dens_plots.ps'
;device, filename=filename, /encapsulated, $
 ; ysize = 9, xsize = 2.5*n_elements(files), /inches

!p.multi = [0,n_elements(files),5,0,1]
!p.charthick = 2.5
!p.charsize = 1.2
!p.thick = 2.5
!x.margin = [14,3]
multiplot

for i = 0, n_elements(files)-1 do begin
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
    ;align, h, g, d, s, 15.
    
;select only the gas particles in the disk
    ind = where(sqrt(g.x^2+g.y^2) lt 15 and abs(g.z) lt 4 and g.tempg lt 15000)
    g = g[ind]

    if(keyword_set(zlimit) ne 0) then begin 
       ind2 = where(s.z lt zlimit and s.z gt -zlimit)
       s = s[ind2]
       rfinal = sqrt(s.x^2 + s.y^2)
       angmom, s, jvec, lvec, jtot, ltot
       rdfloat, files[i]+'.rotcur', r_rot, vc_total
       sort = sort(rfinal)
       rfinal = rfinal[sort]
       sort = sort(sort)
       jzmax = spline(r_rot, vc_total, rfinal)*rfinal
       jzmax = jzmax[sort]
       rfinal = rfinal[sort]
       ind = where(jvec[*,2]/jzmax gt 0.6 and jvec[*,2]/jzmax lt 1.2)
       s = s[ind]
    endif

    g_prof = prof(g, 'gas', h.time, nbins = nbins, rmax = rmax)
    s_prof = prof(s, 'star', h.time, nbins = nbins, rmax = rmax)

    sfr_thresh = fltarr(nbins)+sfrthresh ; sfr threshold for this run in amu/cc
    
; fitpar# is the array p of initial guesses, with components:
; p[0] = break radius (not required for the fit, but returned for practicality)
; p[1] = sfc brightness of inner expo
; p[2] = scale radius of inner expo
; p[3] = sfc brightness of outer expo
; p[4] = scale radius of outer expo

    fitpar =  [0,0,0,0,0]
    
; ********************************************************
; make the fits - select the region you want to fit first
; ********************************************************

    ind = where(s_prof.rbins gt rmin[i] and s_prof.rbins lt rclip[i])
    
    dblexpfit, s_prof.rbins[ind], s_prof.rho[ind], $
      s_prof.rho[ind]/sqrt(s_prof.num_in_bin[ind]), fitpar, /no_guess

; set the densities and masses to physical units

    s_prof.rho = s_prof.rho * msol              ; Msol/kpc^2
    s_prof.sfr = s_prof.sfr*msol                ; Msol/kpc^2/yr 
    g_prof.norm_rho = SYS_AMUCC*g_prof.norm_rho ; amu/cm^3
    g_prof.rho = g_prof.rho*msol*1e-6           ; Msol/pc^2
; **************** 
; STELLAR DENSITY
; ****************

    if(keyword_set(sfrtitle)) then $
      title = textoidl('\rho_{thresh} = ')+string(sfrthresh,format='(F6.2)') $
    else if (keyword_set(titles)) then title = titles[i] $
    else title = 't = ' + string(h.time, format = '(F5.2)') + ' Gyr'

    if (i eq 0)  then begin
        plot, s_prof.rbins, s_prof.rho, /ylog, $
          ytitle = textoidl('\Sigma_{star}')+"[M"+sunsymbol()+" kpc!u-2!n]", $
          title = title,$  
          ystyle=1, yrange=[5e4,1e11],$
          xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
    endif else $ 
      plot, s_prof.rbins, s_prof.rho, /ylog, $
      title = title, $
      ystyle=1, yrange=[5e4,1e11],$
      xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5


    r1 = findgen(100)/100.*(fitpar[0])+3.
    r2 = findgen(100)/100.*(rmax-fitpar[0]+3.) + fitpar[0]-3.
    oplot, r1, msol*fitpar[1]*exp(-r1/fitpar[2]), line = 1
    oplot, r2, msol*fitpar[3]*exp(-r2/fitpar[4]), line = 1
    oplot, [fitpar[0],fitpar[0]], [9.9e4,1e11], line = 3

    multiplot

; *******************
; GAS VOLUME DENSITY
; *******************

    if (i eq 0)  then begin
        plot, g_prof.rbins, g_prof.norm_rho, /ylog, $
          ystyle = 1, yrange=[2e-3,5e2],$
          xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1,$
          ytitle = textoidl('\rho_{gas} [amu cm^{-3}]')
    endif else $
      plot, g_prof.rbins, g_prof.norm_rho, /ylog, $
      ystyle = 1, yrange=[2e-3,5e2],$
      xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1
    
    oplot, g_prof.rbins, g_prof.norm_rho + g_prof.stddev_rho, line = 1
    oplot, g_prof.rbins, g_prof.norm_rho - g_prof.stddev_rho, line = 1
    oplot, g_prof.rbins, sfr_thresh, line = 2
    oplot, [fitpar[0],fitpar[0]], [.0001,1e6], line = 3
    
    multiplot


; ********************
; GAS SURFACE DENSITY
; ********************
    if (i eq 0) then begin
       plot, g_prof.rbins, g_prof.rho, /ylog, $
             xticks = 2, xtickinterval = 5, xrange = xrange, xstyle = 1, yrange = [1e0,1e2], $
             ytitle = textoidl('\Sigma_{gas}')+" [M"+sunsymbol()+" pc!u-2!n]"
    endif else $
       plot, g_prof.rbins, g_prof.rho, /ylog, $
       xticks = 2, xtickinterval = 5, xrange = xrange, xstyle = 1, yrange = [1e0,1e2]
    
    oplot,  [fitpar[0],fitpar[0]], [1e-5,1e5], line = 3

    multiplot

; ********************
; STAR FORMATION RATE
; ********************

    if (i eq 0) then begin 
       plot, s_prof.rbins, s_prof.sfr, /ylog, $
             xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1, yrange = [1e-4,5e0], $
             ytitle = textoidl(' \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]"
    endif else $
       plot, s_prof.rbins, s_prof.sfr, /ylog, $
             xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1, yrange = [1e-4,5e0]

    
    oplot, [fitpar[0],fitpar[0]], [1e-5,1e5], line = 3
    
    multiplot

; *********
; MEAN AGE
; *********

    if(i eq n_elements(files)-1) then xrange=[0,rmax]
    
    if (i eq 0) then begin
        plot, s_prof.rbins, (h.time-s_prof.age)/h.time, $
              xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1,$
              xtitle = 'r [kpc]', ystyle = 1, $
              ytitle = 'Mean age [Gyr/t]', yrange = [.2,.95]
    endif else $
      plot, s_prof.rbins, (h.time-s_prof.age)/h.time, $
            xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1, $
            xtitle = 'r [kpc]', ystyle = 1, yrange = [.2,.95]

    oplot, [fitpar[0],fitpar[0]], [0,1], line = 3
    
    multiplot 


endfor


;device, /close

multiplot, /reset
multiplot, /default
!p.multi = 0
set_plot, 'x'

end
