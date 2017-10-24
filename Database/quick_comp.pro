;Here a quick comparison to see how my magnitudes compare with Andrews
pro quick_comp
;msol = 2.362e5
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
loadct,39

mymagfile = '../Fabio_Runs/h1201-h1755-X5X5g1bwK/h1201-h1755-X5X5g1bwK.00084.amiga_vir.halos.star99.Mv.fits'
andysmagfile = 'Data_files/X5mags.siman.sb99'

mymags = MRDFITS(mymagfile,1)
ind = where(mymags.id ne -1)
mymags = mymags[ind]
readcol,andysmagfile,id,b,v,k,format='I,F,F,F'

match = fltarr(2,N_ELEMENTS(id))
ct_match = 0
FOR ct = 0, N_ELEMENTS(id) - 1 DO BEGIN
   index = where(mymags.id eq id[ct])
   IF (index[0] ne -1) THEN BEGIN
       match[*,ct_match] =[index[0],ct] 
       ct_match = ct_match+1
       ;print,"Matching ",id[ct],' to ',mymags[index[0]].id,' (',ct,', ',index[0],')'
   ENDIF; ELSE print,"Matching ",id[ct]," but no luck"
ENDFOR
   match = match[*,0:ct_match-1]

set_plot,'x'
plot,mymags[match[0,*]].v,v[match[1,*]] ,psym = 3,xrange = [-23,-14],yrange=[-23,-14],xtitle = 'My Magnitudes, V',ytitle = 'Andys Magnitudes, V',title = 'Magnitude Comparison, V'
oplot,[ -23,-14],[-23,-14]
set_plot,'ps'
device, filename='../figures/compX5X5.eps'
plot,mymags[match[0,*]].v,v[match[1,*]] ,psym = 3,xrange = [-25,-8],yrange=[-25,-8],xtitle = 'My Magnitudes, V',ytitle = 'Andys Magnitudes, V',title = 'Magnitude Comparison, V'
oplot,[ -25,-8],[-25,-8]
device,/close
stop
plot,mymags[match[0,*]].b - 4 , b[match[1,*]] ,psym = 3,xrange = [-25,-8],yrange=[-25,-8],xtitle = 'My Magnitudes, B',ytitle = 'Andys Magnitudes, B',title = 'Magnitude Comparison, B'
oplot,[ -25,-8],[-25,-8]
stop
plot,mymags[match[0,*]].k - 1.0,k[match[1,*]] ,psym = 3,xrange = [-25,-8],yrange=[-25,-8],xtitle = 'My Magnitudes, K',ytitle = 'Andys Magnitudes, K',title = 'Magnitude Comparison, K'
oplot,[ -25,-8],[-25,-8]
;device,/close
stop
end
