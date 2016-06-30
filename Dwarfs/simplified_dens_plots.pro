; make some star and gas density plots for run 12M
;example: 
;ytitle_array = ['13M','12M','11M','10M','9M']
;title_array = ['100,000','10,000','1000','100','50']
;simp_dens_plots,"dens_fullgrid.txt",25,10,titles=title_array,ytitle=ytitle_array,filename='mass_grid.eps'

;ytitle_array = ['13M','12M','11M','10M','9M']
;title_array = ['100,000','10,000','1000']
;simp_dens_plots,"dens_fullgrid_kroupa.txt",25,10,titles=title_array,ytitle=ytitle_array,filename='../results/massdens_grid_k.eps'


;ytitle_array = ['8kpc','2kpc','1kpc','500pc','200pc','50pc']
;title_array = ['10000','1000']
;title_array = ['100,000','10,000','1000']
;simp_dens_plots,"dens_specgrid_12M_kroupa.txt",50,10,titles=title_array,ytitle=ytitle_array,filename='../results/spec_grid_12M_kroupa.eps'
;simp_dens_plots,"dens_specgrid_10M_kroupa.txt",50,10,titles=title_array,ytitle=ytitle_array,filename='../results/spec_grid_10M_kroupa.eps'



pro simp_dens_plots, list, nbins, rmax, $
                filename = filename, sfrtitle = sfrtitle, titles = titles, ytitles = ytitles, zlimit = zlimit,linestyle = linstyle

msol = 2.362e5
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/cc


readcol, list, format = 'A60,F3,F3', files, rmin, rclip

xrange = [0,14.99]

if (keyword_set(titles))then column = N_ELEMENTS(titles) else column = 5
if (keyword_set(ytitles)) then row = N_ELEMENTS(ytitles) else row = 3

!p.multi = [0,column,row,0,1]
!p.charthick = 2.5
!p.charsize = 1.2
!p.thick = 2.5
!x.margin = [14,3]

set_plot, 'ps'
if(keyword_set(filename) eq 0) then filename = 'dens_plots.ps'
device, filename=filename, /encapsulated,  ysize = 1.5*row, xsize =2.5*column, /inches

;set_plot,'x'

multiplot


for i = 0, n_elements(files)-1 do begin
    print,files[i]
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot
;    align, h, g, d, s, 15.
    
;select only the gas particles in the disk
;    ind = where(sqrt(g.x^2+g.y^2) lt 15 and abs(g.z) lt 5 and g.tempg lt 15000)
;    g = g[ind]
;My galaxies don't necesarily have a disk so I won't select it out.

    if(keyword_set(zlimit) ne 0) then begin 
       ind2 = where(s.z lt zlimit and s.z gt -zlimit)
       s = s[ind2]

       rfinal = sqrt(s.x^2 + s.y^2)
       angmom, s, jvec, lvec, jtot, ltot
       
       rdfloat, files[i]+'.rotcur', r_rot, vc_total
       
       sort = sort(rfinal)
       
       rfinal = rfinal[sort]
       
       jzmax = spline(r_rot, vc_total, rfinal)*rfinal
       
       sort = sort(sort)
       jzmax = jzmax[sort]
       rfinal = rfinal[sort]

       ind = where(jvec[*,2]/jzmax gt 0.6 and jvec[*,2]/jzmax lt 1.2)
       s = s[ind]

    endif


;The prof procedure will generate a profile using cylidrical shells
;perpendicular to the plane of the galaxy.
    g_prof = prof(g, 'gas', h.time, nbins = nbins, rmax = rmax)
    s_prof = prof(s, 'star', h.time, nbins = nbins, rmax = rmax)

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
    print,i

 ;   stop

    if(keyword_set(sfrtitle)) then $
      title = textoidl('\rho_{thresh} = ')+string(sfrthresh[i],format='(F6.2)') $
    else if (keyword_set(titles)) then title = titles[i/row] $
    else title = 't = ' + string(h.time, format = '(F5.2)') + ' Gyr'

    if (keyword_set(ytitles)) then ytitle = ytitles[i MOD row] $
    else ytitle = textoidl('\Sigma_{gas}')+" [M"+sunsymbol()+" pc!u-2!n]"

;     if ((i mod row) eq 0) then begin
;         if (i eq 0) then begin
;             plot, s_prof.rbins, s_prof.rho, /ylog, $
;               ytitle = ytitle, $
;textoidl('\Sigma_{star}')+"[M"+sunsymbol()+" kpc!u-2!n]", $
;             title = title,$  
;               ystyle=1, yrange=[1e4,1e11],$
;               xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
;         endif else begin
;             plot, s_prof.rbins, s_prof.rho, /ylog,title = title,ystyle=1, yrange=[1e4,1e11],xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
;         endelse
;     endif else begin
;         if ((i/row) eq 0) then begin
;             plot, s_prof.rbins, s_prof.rho, /ylog,ystyle=1, yrange=[1e4,1e11], xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5,ytitle = ytitle
;         endif else begin
;             plot, s_prof.rbins, s_prof.rho, /ylog,$              $
;               ystyle=1, yrange=[1e4,1e11],$
;               xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
;         endelse
;     endelse

;     r1 = findgen(100)/100.*(fitpar[0])+3.
;     r2 = findgen(100)/100.*(rmax-fitpar[0]+3.) + fitpar[0]-3.
;     oplot, r1, msol*fitpar[1]*exp(-r1/fitpar[2]), line = 1
;     oplot, r2, msol*fitpar[3]*exp(-r2/fitpar[4]), line = 1
;     oplot, [fitpar[0],fitpar[0]], [1e4,1e11], line = 3

;     multiplot

; *******************
; GAS VOLUME DENSITY
; *******************

 ;   if (i eq 0)  then begin
 ;       plot, g_prof.rbins, g_prof.norm_rho, /ylog, $
 ;         ystyle = 1, yrange=[2e-3,5e2],$
 ;         xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1,$
 ;         ytitle = textoidl('\rho_{gas} [amu cm^{-3}]')
 ;   endif else $
 ;     plot, g_prof.rbins, g_prof.norm_rho, /ylog, $
 ;     ystyle = 1, yrange=[2e-3,5e2],$
 ;     xticks = 2, xtickinterval = 5, xrange=xrange, xstyle=1
 ;   
 ;   oplot, g_prof.rbins, g_prof.norm_rho + g_prof.stddev_rho, line = 1
 ;   oplot, g_prof.rbins, g_prof.norm_rho - g_prof.stddev_rho, line = 1
 ;   oplot, [fitpar[0],fitpar[0]], [.0001,1e6], line = 3
 ;   
 ;   multiplot


; ********************
; GAS SURFACE DENSITY
; ********************
    if ( (i mod row) eq 0) then begin
        if (i eq 0) then begin
        plot, g_prof.rbins, g_prof.rho, /ylog, $
          xticks = 2, xtickinterval = 5, xrange = xrange, xstyle = 1, yrange = [1e-1,1.5e2], $
          ytitle = ytitle,title = title
        endif else begin
        plot, g_prof.rbins, g_prof.rho, /ylog, $
             xticks = 2, xtickinterval = 5, xrange = xrange, xstyle = 1, yrange = [1e-1,1.5e2],title = title
        endelse
    endif else begin
        if ((i/row) eq 0) then begin
        plot, g_prof.rbins, g_prof.rho, /ylog, $
      xticks = 2, xtickinterval = 5, xrange = xrange, xstyle = 1, yrange = [1e-1,1.5e2], $
      ytitle = ytitle

        endif else begin
          plot, g_prof.rbins, g_prof.rho, /ylog, $
       xticks = 2, xtickinterval = 5, xrange = xrange, xstyle = 1, yrange = [1e-1,1.5e2]
       endelse   
    endelse
    oplot,  [fitpar[0],fitpar[0]], [1e-5,1e5], line = 3
    multiplot
 endfor
;stop

;device, /close

multiplot, /reset
multiplot, /default
!p.multi = 0
;set_plot, 'x'

end

