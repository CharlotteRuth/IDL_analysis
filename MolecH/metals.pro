;
;
;
;
;
;
;
;
;
;   Make some metalicity plots of stars & gas
;
;   parts taken from Chris Brook's halo.pro
;
;
;
;
;metals,'/astro/net/nbody1/christensen/MolecH/MWHR/12M_hr/12M_hr.01000.dir/12M_hr.01000',1,2.362e5, /do_plots
;metals,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK.BUG/h603.cosmo50cmb.2304g2bwK.00512/h603.cosmo50cmb.2304g2bwK.00512.halo.1',50000.,1.84793e16,/do_plots
;metals,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG/h603.cosmo50cmb.2304g5bwK.00512.halo.1',50000.,1.84793e16,/do_plots,/ascii
;metals,'/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG/h603.cosmo50cmb.2304g5bwK.00512.halo.1',50000.,1.84793e16,/do_plots,/ascii

pro metals, filename, dKpcUnit, dMsolUnit, fehs, fehg, ofes, ofeg, gas, h=h, g=g,s=s,$
            DISK_ONLY = disk_only, DO_PLOTS = do_plots,ascii = ascii

if(keyword_set(h) eq 0) then begin
    rtipsy, filename, h, g, d, s
; center the system
    cofm, h, g, d, s, /pot
endif


; define some constants
ms_units=2.325e5
v_units=1
dist_units=1
t_units=SQRT((dKpcUnit*3.086d21)^3/(6.67d-8*dMsolUnit*1.99d33))/(3600.*24.*365.24)
XSOLfe=0.117E-2
XSOLO=0.96E-2
XSOLH=0.706

; read in the data files
oxfile = filename+'.OxMassFrac'
fefile = filename+'.FeMassFrac'
HIfile = filename+'.HI'
HeIfile = filename+'.HeI'
HeIIfile = filename+'.HeII'


readarr, oxfile, h, ox, ascii = ascii
readarr, fefile, h, fe, ascii = ascii
readarr, HIfile, h, HI, ascii = ascii
readarr, HeIfile, h, HeI, ascii = ascii
readarr, HeIIfile, h, HeII, ascii = ascii

;readcol, oxfile, ox
;readcol, fefile,  fe
;readcol, HIfile,  HI
;readcol, HeIfile,  HeI
;readcol, HeIIfile,  HeII
;ox = ox[1:N_ELEMENTS(ox) -1]
;fe = fe[1:N_ELEMENTS(ox) -1]
;HI = HI[1:N_ELEMENTS(ox) -1]
;HeI = HeI[1:N_ELEMENTS(ox) -1]
;HeII = HeII[1:N_ELEMENTS(ox) -1]

; make new arrays of structures that contain metalicity data

gas = replicate({mass: 1.,x: 1.,y : 1., z:1.,$
                 vx:1.,vy:1.,vz:1.,$
                 dens:1.,tempg:1., h: 1. , zmetal: 1., phi: 1.,$
                ox: 1., fe: 1., HI: 1., HeI: 1., HeII: 1.},h.ngas)

gas_ind = findgen(h.ngas)
dark_ind = findgen(h.ndark)+h.ngas
stars_ind = findgen(h.nstar)+h.ngas+h.ndark

; set up the gas
gas.x =       g.x*dKpcUnit
gas.y =       g.y*dKpcUnit 
gas.z =       g.z*dKpcUnit
gas.vx =      g.vx
gas.vy =      g.vy 
gas.vz =      g.vz
gas.mass =    g.mass
gas.dens =    g.dens
gas.tempg =   g.tempg
gas.h =       g.h
gas.zmetal =  g.zmetal
gas.phi =     g.phi
gas.ox =      ox[gas_ind]
gas.fe =      fe[gas_ind]
gas.HI =      HI[gas_ind]
gas.HeI =    HeI[gas_ind]
gas.HeII =  HeII[gas_ind]

; free the array to conserve memory
;g = 0

; set up the stars
stars = replicate({mass: 1.,x: 1.,y : 1., z:1.,$
                   vx:1.,vy:1.,vz:1.,metals:1.,tform:1.,eps: 1.,phi: 1.,$
                  ox: 1., fe: 1., HI: 1., HeI: 1., HeII: 1.},h.nstar)
stars.x =      s.x*dKpcUnit
stars.y =      s.y*dKpcUnit 
stars.z =      s.z*dKpcUnit
stars.vx =     s.vx
stars.vy =     s.vy 
stars.vz =     s.vz
stars.mass =   s.mass
stars.metals = s.metals
stars.tform =  s.tform
stars.eps =    s.eps
stars.phi =    s.phi
stars.ox =     ox[stars_ind]
stars.fe =     fe[stars_ind]
stars.HI =     HI[stars_ind]
stars.HeI =   HeI[stars_ind]
stars.HeII = HeII[stars_ind]

; free the array
;s = 0
;d = 0 ; won't need the dark matter

;
; RADIAL POSITIONS (ASSUMING X-Y IS THE DISK PLANE)
;

rs = sqrt(stars.x^2+stars.y^2)
rg = sqrt(gas.x^2 + gas.y^2)

;
; ONLY KEEP THE GAS PARTICLES THAT ARE IN THE DISK? 
;

if(keyword_set(DISK_ONLY)) then begin
    gaskeep = where(rg lt disk_only and abs(gas.z) lt 1.0)
    gas = gas[gaskeep]
    rg = rg[gaskeep]
endif

; some stars/gas particles have zero metalicity, so set their metallicity to lowest
; non-zero metalicity

mets = where(stars.metals gt 0.,comp=metsfree)
metg = where(gas.zmetal gt 0., comp=metgfree) 

fehs = fltarr(h.nstar)
oxhs = fltarr(h.nstar)
ofes = fltarr(h.nstar)
fehg = fltarr(n_elements(gas))
oxhg = fltarr(n_elements(gas))
ofeg = fltarr(n_elements(gas))

; **************
; Iron fraction 
; **************
    
fehs[mets] = stars[mets].fe/$
  (1.-stars[mets].metals-(stars[mets].HeI+stars[mets].HeII))
fehs[metsfree] = min(fehs[mets])

lowmet = where(fehs lt 1e-5*XSOLFe/XSOLH, count)

print, '# with low metalicity =', count

fehs[lowmet] = 1e-5*XSOLFe/XSOLH

fehg[metg] = gas[metg].fe/$
  (1.-gas[metg].zmetal-(gas[metg].HeI+gas[metg].HeII))

; **************
; Iron fraction 
; **************
    
oxhs[mets] = stars[mets].ox/$
  (1.-stars[mets].metals-(stars[mets].HeI+stars[mets].HeII))
oxhs[metsfree] = min(oxhs[mets])

lowmet = where(oxhs lt 1e-5*XSOLO/XSOLH, count)

print, '# with low metalicity =', count

oxhs[lowmet] = 1e-5*XSOLO/XSOLH

oxhg[metg] = gas[metg].ox/$
  (1.-gas[metg].zmetal-(gas[metg].HeI+gas[metg].HeII))


; ******************
; Oxygen/Iron ratio
; ******************

ofes[mets] = stars[mets].ox/stars[mets].fe
ofes[metsfree] = min(ofes[mets])

ofeg[metg] = gas[metg].ox/gas[metg].fe

;fehs[metsfree] = min(fehs[mets])
print, alog10(min(fehs[mets]))-alog10(XSOLfe/XSOLH)
if(metgfree[0] ne -1) then fehg[metgfree] = min(fehg[metg])

ofes[metsfree] = min(ofes[mets])
if(metgfree[0] ne -1) then ofeg[metgfree] = min(ofeg[metg])


if(keyword_set(do_plots)) then begin
    set_plot,'x'
    set_plot, 'ps'
    
    loadct, 39

    device, filename = filename+'.metals.ps', /color, ysize = 5, xsize = 9, /inches
    !p.charsize = 1.2
    !p.charthick = 2.0
    !p.thick = 2.0
;
; METALICITY HISTOGRAMS
;

    !p.multi = [0,2,1]

    age = (MAX(s.tform) - s.tform)*t_units/1e9
    
    age1 = where(age lt 0.1)
    age2 = where(age gt 1.0 and age lt 3.0)
    age3 = where(age gt 3.0)

    binsize = 0.5

    stddevg = 0.0

    masshistg = hist1d(rg, gas.mass, binsize = binsize)
    meanfeg   = hist1d(rg, fehg*gas.mass, binsize = binsize, $
                       obin = xg, stddev = stddevg)
    meanfeg   = meanfeg/masshistg
    meanfehg  = alog10(meanfeg)-alog10(XSOLFe/XSOLH)   

    masshists = hist1d(rs, stars.mass, binsize = binsize)
    meanfes   = hist1d(rs, fehs*stars.mass, binsize = binsize, $
                       obin = xs, stddev = stddev_all)
    meanfes   = meanfes/masshists
    meanfehs  = alog10(meanfes)-alog10(XSOLFe/XSOLH)   

    masshists1 = hist1d(rs[age1], stars[age1].mass, binsize = binsize)
    meanfes1   = hist1d(rs[age1], fehs[age1]*stars[age1].mass, $
                        stddev = stddev1, binsize = binsize, obin = xs1)
    meanfes1   = meanfes1/masshists1
    meanfehs1  = alog10(meanfes1)-alog10(XSOLFe/XSOLH)   

    masshists2 = hist1d(rs[age2], stars[age2].mass, binsize = binsize)
    meanfes2   = hist1d(rs[age2], fehs[age2]*stars[age2].mass, $
                        stddev = stddev2, binsize = binsize, obin = xs2)
    meanfes2   = meanfes2/masshists2
    meanfehs2  = alog10(meanfes2)-alog10(XSOLFe/XSOLH)   

    masshists3 = hist1d(rs[age3], stars[age3].mass, binsize = binsize)
    meanfes3   = hist1d(rs[age3], fehs[age3]*stars[age3].mass, $
                        stddev = stddev3, binsize = binsize, obin = xs3)
    meanfes3   = meanfes3/masshists3
    meanfehs3  = alog10(meanfes3)-alog10(XSOLFe/XSOLH)   
    
    
    meanoxg   = hist1d(rg, oxhg*gas.mass, binsize = binsize, $
                       obin = xg, stddev = stddevg)
    meanoxg   = meanoxg/masshistg
    meanoxhg  = alog10(meanoxg)-alog10(XSOLO/XSOLH)   

    meanoxs   = hist1d(rs, oxhs*stars.mass, binsize = binsize, $
                       obin = xs, stddev = stddev_all)
    meanoxs   = meanoxs/masshists
    meanoxhs  = alog10(meanoxs)-alog10(XSOLO/XSOLH)    

    masshists1 = hist1d(rs[age1], stars[age1].mass, binsize = binsize)
    meanoxs1   = hist1d(rs[age1], oxhs[age1]*stars[age1].mass, $
                        stddev = stddev1, binsize = binsize, obin = xs1)
    meanoxs1   = meanoxs1/masshists1
    meanoxhs1  = alog10(meanoxs1)-alog10(XSOLO/XSOLH)   

    masshists2 = hist1d(rs[age2], stars[age2].mass, binsize = binsize)
    meanoxs2   = hist1d(rs[age2], oxhs[age2]*stars[age2].mass, $
                        stddev = stddev2, binsize = binsize, obin = xs2)
    meanoxs2   = meanoxs2/masshists2
    meanoxhs2  = alog10(meanoxs2)-alog10(XSOLO/XSOLH)   

    masshists3 = hist1d(rs[age3], stars[age3].mass, binsize = binsize)
    meanoxs3   = hist1d(rs[age3], oxhs[age3]*stars[age3].mass, $
                        stddev = stddev3, binsize = binsize, obin = xs3)
    meanoxs3   = meanoxs3/masshists3
    meanoxhs3  = alog10(meanoxs3)-alog10(XSOLO/XSOLH)  

    meanofes = hist1d(rs, ofes*stars.mass, binsize = 0.2)
    meanofes = meanofes/masshists
    meanofes = alog10(meanofes)-alog10(XSOLO/XSOLFe)
           

    plot, xg, meanfehg, $
      xtitle = 'R [kpc]', ytitle = '[Fe/H]', yrange = [-1.5,0.5], ystyle = 1, psym = 10,xrange = [0,20]
    oplot, xs, meanfehs, line = 5
    oplot, xs1, meanfehs1, color = 50, psym = 10
    oplot, xs2, meanfehs2, color = 150, psym = 10
    oplot, xs3, meanfehs3, color = 250, psym = 10
;    legend, /right, ['gas', 'all stars'], $
;      line = [0,5]
    legend, /right, ['gas', 'all stars', 'age < 0.1 Gyr', '1 < age < 3 Gyr', 'age > 5 Gyr'], $
      line = [0,5,0,0,0], color = [0,0,50,150,250]

;    plot, xg, alog10(stddevg) - alog10(XSOLfe/XSOLH), $
;      xtitle = 'R [kpc]', ytitle = textoidl('\sigma_{[Fe/H]}'), psym = 10,xrange = [0,20]
;    oplot, xs, alog10(stddev_all) - alog10(XSOLfe/XSOLH), line = 5
;    oplot, xs1, alog10(stddev1) - alog10(XSOLfe/XSOLH), color = 50, psym = 10
;    oplot, xs2, alog10(stddev2) - alog10(XSOLfe/XSOLH), color = 150, psym = 10
;    oplot, xs3, alog10(stddev3) - alog10(XSOLfe/XSOLH), color = 250, psym = 10
;stop
    plot, xg, meanoxhg, $
      xtitle = 'R [kpc]', ytitle = '[Ox/H]', yrange = [-1.5,0.5], ystyle = 1, psym = 10,xrange = [0,20]
    oplot, xs, meanoxhs, line = 5
    oplot, xs1, meanoxhs1, color = 50, psym = 10
    oplot, xs2, meanoxhs2, color = 150, psym = 10
    oplot, xs3, meanoxhs3, color = 250, psym = 10

    legend, /right, ['gas', 'all stars', 'age < 0.1 Gyr', '1 < age < 3 Gyr', 'age > 5 Gyr'], $
      line = [0,5,0,0,0], color = [0,0,50,150,250]

 ;   plot, xg, alog10(stddevg) - alog10(XSOLO/XSOLH), $
 ;     xtitle = 'R [kpc]', ytitle = textoidl('\sigma_{[Fe/H]}'), psym = 10,xrange = [0,20]
 ;   oplot, xs, alog10(stddev_all) - alog10(XSOLO/XSOLH), line = 5

    device, /close
     set_plot, 'x'
    !p.multi = 0
endif


end
