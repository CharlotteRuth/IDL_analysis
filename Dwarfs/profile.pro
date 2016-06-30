; make some star and gas density plots for run 12M
;example: 
;ytitle_array = ['13M','12M','11M','10M','9M']
;title_array = ['100,000','10,000','1000','100','50']
;simp_dens_plots,"dens_fullgrid.txt",25,10,titles=title_array,ytitle=ytitle_array,filename='mass_grid.eps'

;ytitle_array = ['8kpc','2kpc','1kpc','500pc','200pc','50pc']
;title_array = ['10000','1000']
;simp_dens_plots,"dens_specgrid.txt",50,10,titles=title_array,ytitle=ytitle_array,filename='spec_grid.eps'

PRO deVauc,x,A,F,pder
;x is the radius, 
;A[0]= effective radius, 
;A[1] = surface brightness at the effective radius
F = A[1]*EXP(-7.67*((X/A[0])^0.25 - 1))
d1 = EXP(-7.67*((X/A[0])^0.25 - 1))
d2 = A[1]*EXP(-7.67*((X/A[0])^0.25 - 1))*(-1.9175)*x^0.25*A[0]^(-1.25)
pder =[[d1],[d2]]
END


PRO exp,x,A,F,pder
;x is the radius, 
;A[0]= effective radius, 
;A[1] = surface brightness at the effective radius
F = A[1]*EXP((-1.0)*X/A[0])
d1 = EXP((-1.0)*X/A[0])
d2 = A[1]*EXP(x/A[0])*x*A[0]^(-2.0)
pder =[[d1],[d2]]
END

pro profile, list, nbins, rmax, $
                filename = filename, sfrtitle = sfrtitle, titles = titles, ytitles = ytitles, zlimit = zlimit

msol = 2.362e5
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c

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

if(keyword_set(filename) eq 0) then filename = 'profile.ps'

device, filename=filename, /encapsulated, $
  ysize = 1.5*row, xsize = 2.5*column, /inches 
;set_plot,'x'

multiplot

close,2
openw,2,'fits.dat'

for i = 0, n_elements(files)-1 do begin
    print,files[i]
    rtipsy, files[i], h, g, d, s
    cofm, h, g, d, s, /pot

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
      s_prof.rho[ind]/sqrt(s_prof.num_in_bin[ind]), fitpar, red_chisq=exp2chi, /no_guess
    
;Make a guess at the de Vaucouleurs profile
    A=[1.0,1.0]
    expchi=1e11
    weights = FLTARR(N_ELEMENTS(ind))+1.0
    deVaucFit=CURVEFIT(s_prof.rbins[ind],s_prof.rho[ind],weights,A,sigma,FUNCTION_NAME='exp',CHISQ=expchi)
    expchi=expchi/(n_elements(ind) - n_elements(A) - 1.)
    print,'EXP: ',A,expchi

;Make a guess at the de Vaucouleurs profile
    A=[1.0,1.0]
    dVchi=1e11
    weights = FLTARR(N_ELEMENTS(ind))+1.0
    deVaucFit=CURVEFIT(s_prof.rbins[ind],s_prof.rho[ind],weights,A,sigma,FUNCTION_NAME='deVauc',CHISQ=dVchi)
    dVchi=dVchi/(n_elements(ind) - n_elements(A) - 1.)
    print,'DeVauc: ',A,dVchi

    if(dVchi lt expchi AND dVchi lt exp2chi) then best = '   deVauc' $
    else if (exp2chi lt expchi AND exp2chi lt dVchi) then best = '   2 Exp'$
    else best='   1 Exp'

    printf,2,files[i],exp2chi,expchi,dVchi,best

; set the densities and masses to physical units

    s_prof.rho = s_prof.rho * msol              ; Msol/kpc^2
    s_prof.sfr = s_prof.sfr*msol                ; Msol/kpc^2/yr 
    g_prof.norm_rho = SYS_AMUCC*g_prof.norm_rho ; amu/cm^3
    g_prof.rho = g_prof.rho*msol*1e-6           ; Msol/pc^2
; **************** 
; STELLAR DENSITY
; ****************
    if(keyword_set(sfrtitle)) then $
      title = textoidl('\rho_{thresh} = ')+string(sfrthresh[i],format='(F6.2)') $
    else if (keyword_set(titles)) then title = titles[i/row] $
    else title = 't = ' + string(h.time, format = '(F5.2)') + ' Gyr'

    if (keyword_set(ytitles)) then ytitle = ytitles[i MOD row] $
    else ytitle = textoidl('\Sigma_{star}')+" [M"+sunsymbol()+" pc!u-2!n]"

     if ((i mod row) eq 0) then begin
         if (i eq 0) then begin
             plot, s_prof.rbins, s_prof.rho, /ylog, $
               ytitle = ytitle, title = title,$  
               ystyle=1, yrange=[1e4,1e11],$
               xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
         endif else begin
             plot, s_prof.rbins, s_prof.rho, /ylog,title = title,ystyle=1, yrange=[1e4,1e11],xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
         endelse
     endif else begin
         if ((i/row) eq 0) then begin
             plot, s_prof.rbins, s_prof.rho, /ylog,ystyle=1, yrange=[1e4,1e11], xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5,ytitle = ytitle
         endif else begin
             plot, s_prof.rbins, s_prof.rho, /ylog,$              $
               ystyle=1, yrange=[1e4,1e11],$
               xstyle=1, xrange=xrange, xticks = 2, xtickinterval = 5
         endelse
     endelse

     r1 = findgen(100)/100.*(fitpar[0])+3.
     r2 = findgen(100)/100.*(rmax-fitpar[0]+3.) + fitpar[0]-3.
     oplot, r1, msol*fitpar[1]*exp(-r1/fitpar[2]), line = 1
     oplot, r2, msol*fitpar[3]*exp(-r2/fitpar[4]), line = 1
     oplot, [fitpar[0],fitpar[0]], [1e4,1e11], line = 3

     multiplot
 endfor
;stop
close,2
device, /close

multiplot, /reset
multiplot, /default
!p.multi = 0
set_plot, 'x'

end
