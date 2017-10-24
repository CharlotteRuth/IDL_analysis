; Find V200 for a list of galaxy files
;v200,'halo_list.dat','X2X2tf'

pro mag_test, list, outputfile
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
vcirc2 = fltarr(n_elements(files))
rcirc2 = fltarr(n_elements(files))
v200 = fltarr(n_elements(files))
stellarmass = fltarr(n_elements(files))
baryonicmass = fltarr(n_elements(files))
mag_file = '../Fabio_Runs/h1201-h1755-X2X2g1bwK/h1201-h1755-X2X2g1bwK.00512.amiga.halos.Mv.fits'
halo_mags = MRDFITS(mag_file,1)
num_halos = N_ELEMENTS(halo_mags)

indices = STRTRIM(sindgen(num_halos),2)
imag = fltarr(n_elements(files))
plotted_i = fltarr(n_elements(files))
kmag = fltarr(n_elements(files))
plotted_k = fltarr(n_elements(files))
nbins = 5000
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
    kmag[ct] = halo_mags[halo_index].i
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
    v200[ct] = SQRT(m200*dMsolUnit*msol*grav/r200/dkpcUnit/kpc)

    rcirc = 20./220.*v200/hubble
    mcirc = TOTAL(mass[where(distances*dkpcUnit LE rcirc[ct])])
    vcirc[ct] = SQRT(mcirc*dMsolUnit*msol*grav/rcirc[ct]/kpc)

    velocity = SQRT(tmass*dMsolUnit*msol*grav/bins/dkpcUnit/kpc)


    IF ((Where (density le 200.0))[0] eq -1) THEN vcirc[ct] = 1.0
    stellarmass[ct] = TOTAL(s[WHERE(distancess*dkpcUnit le rcirc[ct])].mass)*dMsolUnit
    baryonicmass[ct] = TOTAL(s[WHERE(distancess*dkpcUnit le rcirc[ct])].mass)*dMsolUnit + TOTAL(g[WHERE(distancesg*dkpcUnit le rcirc[ct])].mass)*dMsolUnit
    print,'V200: ',v200[ct],' Vcirc: ',vcirc[ct]
    print,'M200: ',m200*dMsolUnit,' Mcirc: ',mcirc*dMsolUnit
    print,'R200: ',r200*dkpcUnit,' Rcirc: ',Rcirc[ct]
    print,'Stellar Mass',stellarmass[ct]
    printf,7,halo_index,rcirc[ct]
;    stop
endfor
close,7

set_plot,'x'
;Estimates the luminosity of the galaxy based on a mass to light ratio
;of 4 Msol/Lsol
plot,3.28-2.5*ALOG10(2.0*stellarmass),kmag,psym = 4
;,xtitle = "3.28 - 2.5 LOG(2*stellar mass)",ytitle = "K magnitude"

set_plot, 'ps'
device, filename=outputfile+'test.eps'
plot,3.28-2.5*ALOG10(2.0*stellarmass),kmag
;,psym = 4,xtitle = "3.28 - 2.5 LOG(2*stellar mass)",ytitle = "K magnitude"
device,/close
stop

set_plot, 'x'

end
