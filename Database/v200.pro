; Find V200 for a list of galaxy files
;v200,'halo_list.dat','X2X2tf'

pro v200, list, outputfile
;msol = 2.362e5
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70

readcol, list, format = 'A60', files
vcirc = fltarr(n_elements(files))
rcirc = fltarr(n_elements(files))
stellarmass = fltarr(n_elements(files))
mag_file = '../Fabio_Runs/h1201-h1755-X2X2g1bwK/h1201-h1755-X2X2g1bwK.00512.amiga.halos.Mv.fits'
halo_mags = MRDFITS(mag_file,1)
num_halos = N_ELEMENTS(halo_mags)
indices = STRTRIM(sindgen(num_halos),2)
imag = fltarr(n_elements(files))
plotted_i = fltarr(n_elements(files))
nbins =300
tmass = fltarr(nbins)
density = fltarr(nbins)
loadct,39

close,7
openw,7,'halo_radii.txt'
for ct = 0, n_elements(files)-1 do begin
    print,files[ct]
    halo_num = (strsplit(files[ct],'.',/extract))[4]
    IF (strpos(halo_num,'0') eq 0) THEN halo_num = STRMID(halo_num,1,STRLEN(halo_num)-1)
    halo_index = where(strcmp(halo_num,indices) eq 1)
    imag[ct] = halo_mags[halo_index].i
;    plotted_i[ct] = halo_index
    rtipsy, files[ct], h, g, d, s
    distancesg = SQRT(g.x^2.0 + g.y^2.0 + g.z^2.0)
    distancess = SQRT(s.x^2.0 + s.y^2.0 + s.z^2.0)
    distancesd = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)
    distances = [distancesg,distancess,distancesd]
    massg = g.mass
    masss = s.mass
    massd = d.mass
    mass = [massg,masss,massd]
    maxd = MAX(distances)
    mind = MIN(distances)
    bins = findgen(nbins)*(maxd-mind)/nbins+mind
    FOR i = 0, nbins-1 DO  tmass[i] = TOTAL(mass[where(distances LE bins[i])])
    density = tmass/(4.0/3.0*!PI*bins^3.0)
    temp = MIN(ABS(density  - 200.0),i200)
    r200 = SPLINE(-1.*density[i200-2:i200+2],bins[i200-2:i200+2],-200.0)
    m200 = TOTAL(mass[where(distances LE r200)])
;    temp = MIN(ABS(density  - 200.0),i200)
;    r200 = bins[i200]
;    m200 = tmass[i200]
;    plot,bins*50000,density,/ylog,yrange=[100,1e9],title = 'Baryonic Mass Profile',xtitle='Radius (kpc)',ytitle = 'Total Mass Density (* p_crit)',/xlog
;    oplot,[MIN(bins)*50000,MAX(BINS)*50000],[200,200],linestyle = 1
    v200 = SQRT(m200*dMsolUnit*msol*grav/r200/dkpcUnit/kpc)
;    stellarmass[ct] = TOTAL(s[WHERE(distancess le r200)].mass)*dMsolUnit 
;    IF ((Where (density le 200.0))[0] eq -1) THEN V200[ct] = 1.0   

    rcirc[ct] = 20./220.*v200/hubble
    mcirc = TOTAL(mass[where(distances*dkpcUnit LE rcirc[ct])])
    vcirc[ct] = SQRT(mcirc*dMsolUnit*msol*grav/rcirc[ct]/kpc)

    IF ((Where (density le 200.0))[0] eq -1) THEN vcirc[ct] = 1.0
    stellarmass[ct] = TOTAL(s[WHERE(distancess*dkpcUnit le rcirc[ct])].mass)*dMsolUnit
    print,'V200: ',v200,' Vcirc: ',vcirc[ct]
    print,'M200: ',m200*dMsolUnit,' Mcirc: ',mcirc*dMsolUnit
    print,'R200: ',r200*dkpcUnit,' Rcirc: ',Rcirc[ct]
    print,'Stellar Mass',stellarmass[ct]
    printf,7,halo_index,rcirc[ct]
endfor
close,7

set_plot,'x'
plot,stellarmass,ALOG10(vcirc),title = 'Simple TF',ytitle = 'Log(V200)', xtitle = 'Stellar Mass',psym = 4,/ylog,/xlog
set_plot, 'ps'
device, filename=outputfile+'stellar.eps'
plot,stellarmass,ALOG10(vcirc),title = 'Simple TF',ytitle = 'Log(V200)', xtitle = 'Stellar Mass',psym = 4,/ylog,/xlog
device,/close


!Y.style = 1
!X.style = 1
mags = findgen((24.0-15.0)/0.1)*0.1 -24.0
slope = -0.128412
velocities = slope*(mags + 15.5) + 1.5
;Fit given in Eke, Navarro and Steimetz 2001

set_plot,'x'
plot,imag - 5*ALOG10(hubble),ALOG10(vcirc),title = 'TF',xtitle = 'MI - 5Log(h)', ytitle = 'Log(V200)',psym = 4,yrange=[1.5,2.8],xrange = [-15,-24]
oplot,mags,velocities
set_plot, 'ps'
device, filename=outputfile+'.eps'
plot,imag - 5*ALOG10(hubble),ALOG10(vcirc),title = 'TF',xtitle = 'MI - 5Log(h)', ytitle = 'Log(V200)',psym = 4,yrange=[1.5,2.8],xrange = [-15,-24]
oplot,mags,velocities
device,/close

set_plot, 'x'

end
