pro rotCurve_tipsy_master,file,output_graph = output_graph, nbins = nbins

;NAME:
;  rotCurve
;
;PURPOSE:
;   This function will plot the rotation curve for a tipsy file
;
;
;CALLING SEQUENCE:
;   IDL>rotCurve,"myfile.std"
;
;REQUIRED INPUTS:
;   "myfile.std" = Name of a TIPSY snapshot file
;
;OPTIONAL INPUTS:
;
;   output_graph = If you would like the rotation curve and the
;   density profile saved to a file, the is the root of that file
;
;   nbins = is the number of density bins
;

;OUTPUTS:
;  A  graph of the rotation curve
;
;EXTERNAL ROUTINES:
;   rtipsy

;EXTERNAL DATA:
;
;COMMENTS:
;
;REVISION HISTORY:
;

loadct,39
msol = 2d33 ; mass of Sum in grams /1e16
dMsolUnit = 2.362e5 ; Solar mass in system units
kpc = 3.085d21 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 1.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
IF (keyword_set(nbins)) THEN nbins = nbins ELSE nbins = 100
IF (keyword_set(maxdistance)) THEN maxdistance = maxdistance ELSE maxdistance = 1000

bind = maxdistance/dKpcUnit/nbins
bins = (findgen(nbins)+1)*bind
tmass = fltarr(nbins)
density = fltarr(nbins)

base = '/astro/net/scratch1/christensen/DwarfResearch/ResMassTests/'
files10 = base + ['1E6R/12M_k/o12M_1.00300','1E5R/12M_k/o12M_1.00300','1E4R/12M_k/o12M_1.00300','1E3R/12M_k/o12M_1.00300','1E2R/12M_k/o12M_1.00300','5E1R/12M_k/o12M_1.00300']
for ct = 0, N_ELEMENTS(files10) - 1 DO rotCurve,files10[ct],iteration = ct
legend, ['1e6','1e5','1e4','1000','100','50'],color = colors,linestyle = [0,0,0,0,0,0] ,/right
stop

files12 = base + ['1E6R/10M_k/o10M_1.00300','1E5R/10M_k/o10M_1.00300','1E4R/10M_k/o10M_1.00300','1E3R/10M_k/o10M_1.00300','1E2R/10M_k/o10M_1.00300','5E1R/10M_k/o10M_1.00300']
for ct = 0, N_ELEMENTS(files10) - 1 DO rotCurve,files12[ct],iteration = ct
legend, ['1e6','1e5','1e4','1000','100','50'],color = colors,linestyle = [0,0,0,0,0,0] ,/right
END

pro rotCurve_tipsy,file,iteration = iteration
loadct,39
msol = 2d33 ; mass of Sum in grams /1e16
dMsolUnit = 2.362e5 ; Solar mass in system units
kpc = 3.085d21 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 1.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
IF (keyword_set(nbins)) THEN nbins = nbins ELSE nbins = 100
IF (keyword_set(maxdistance)) THEN maxdistance = maxdistance ELSE maxdistance = 1000

bind = maxdistance/dKpcUnit/nbins
bins = (findgen(nbins)+1)*bind
tmass = fltarr(nbins)
density = fltarr(nbins)

colors=[240,190,150,100,60,20]
rtipsy, file, h, g, d, s
distancesg = SQRT(g.x^2.0 + g.y^2.0 + g.z^2.0)
IF N_ELEMENTS(s) gt 0 THEN distancess = SQRT(s.x^2.0 + s.y^2.0 + s.z^2.0)
 distancesd = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)
IF N_ELEMENTS(s) gt 0 THEN distances = [distancesg,distancess,distancesd] ELSE distances = [distancesg,distancesd]
massg = g.mass*dMsolUnit
IF N_ELEMENTS(s) gt 0 THEN masss = s.mass*dMsolUnit
massd = d.mass*dMsolUnit
IF N_ELEMENTS(s) gt 0 THEN mass = [massg,masss,massd] ELSE mass = [massg,massd]
maxd = MAX(distances)
mind = MIN(distances)

FOR i = 0, nbins-1 DO BEGIN
    ind = where(distances LE bins[i])
    if (ind[0] ne -1)then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
ENDFOR

density = tmass/(4.0/3.0*!PI*bins^3.0)
temp = MIN(ABS(density  - 200.0),i200)
IF(i200 LT 2) THEN i200 = 2
IF(i200 GT nbins -2) THEN i200 = nbins - 2
print,file

velocitycurve = SQRT(tmass*msol*grav/bins/dkpcUnit/kpc)
set_plot,'x'
IF iteration eq 0 Then plot,bins*dKpcUnit,velocitycurve,xtitle = 'Radius',ytitle = 'Velocity (km/s)',title = 'Rotation Curve '+file 
oplot,bins*dKpcUnit,velocitycurve,color = colors[iteration]

IF (keyword_set(output_graph)) THEN BEGIN
    set_plot,'ps'
    device,filename = output_graph+'_RC.pro'
    plot,bins*dKpcUnit,velocitycurve,xtitle = 'Radius',ytitle = 'Velocity (km/s)',title = 'Rotation Curve '+file
    set_plot,'x'
ENDIF
END
