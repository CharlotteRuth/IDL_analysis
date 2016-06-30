PRO halo_info,filenames,fileroot,haloid,validids,munit,lunit,debug = debug
;print,filenames,fileroot,haloid,validids,munit,lunit
loadct,39
formatplot,outplot = outplot
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
H_per_gm = 5.9753790e+23
grav = 6.67e-8
zsolar  =  0.0130215
dtime = 1e8 ;100 MYR
rhounit=munit * gm_per_msol * H_per_gm/lunit^3/cm_per_kpc^3
tunit=sqrt((lunit*3.086d21)^3/(6.67d-8*munit*1.99d33))/(3600.*24.*365.24)

IF keyword_set(outplot) THEN BEGIN
    fgcolor = 0 
    bgcolor = 255
    xsize = 6.25
    ysize = 10
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    xsize = 500
    ysize = 800
ENDELSE

openw,1,fileroot + '.' + haloid+ '_prof.txt'
printf,1,'Filename ','Halo ID','Redshift ','V_circ [km/s]','Virial Mass [Msun]','Dark Matter Mass [Msun]','Stellar Mass [Msun]','Gas Mass [Msun]','HI Mass [Msun]','H2 Mass [Msun]','Cold Gas Mass [Msun]','Hot Gas Mass [Msun]','Gas Metals [Msun]','Cold Gas Metals [Msun]','Hot Gas Metals [Msun]','Half Gas Mass Rad [kpc]','Half Stellar Mass Rad [kpc]',format='(A40,16A28)'
close,1
;IF NOT keyword_set(filebase) THEN filebase = filenames
FOR i = 0, n_elements(validids) - 1 DO BEGIN
    print,filenames[i],', ',validids[i]
 ;   rtipsy,filenames[i],h,g,d,s
;    strel = strsplit(validids[i],'.',/EXTRACT)
;    halo = strel[n_elements(strel) - 1]
    strel = strsplit(filenames[i],'/',/EXTRACT)
    filebase = strel[n_elements(strel) - 1]
    stat = read_stat_struc_amiga(filenames[i] + '.amiga.stat')
    statind = where(stat.group EQ validids[i])
    tipsysatshi,filenames[i],validids[i],lunit,munit,/cutout_rad,num = h,dark = d,gas = g,stars = s,ind = ind,/notipsy
    rtipsy,filenames[i],hAll,gAll,dAll,sAll,/justhead
    read_tipsy_arr,filenames[i] + '.' +'HI',hAll,HI,type = 'float'
    read_tipsy_arr,filenames[i] + '.' +'H2',hAll,H2,type = 'float'
    HI = HI[ind]
    H2 = H2[ind]
    HI = HI[0:h.ngas - 1]
    H2 = H2[0:h.ngas - 1]

    gmass_HI = total(g.mass*HI)*munit
    gmass_H2 = total(g.mass*H2)*munit
    cold = where(g.tempg LT 2e4 AND g.dens*rhounit GT 0.01, complement = hot)
    gmass_cold = total(g[cold].mass)*munit
    gmass_coldz = total(g[cold].mass*g[cold].zmetal)*munit
    gmass_hot = total(g[hot].mass)*munit
    gmass_hotz = total(g[hot].mass*g[hot].zmetal)*munit
    gmassz = total(g.mass*g.zmetal)*munit
;    read_tipsy_arr,filenames[i] + '.FeMassFrac',h,fe, part = 'gas'
;    read_tipsy_arr,filenames[i] + '.OxMassFrac',h,ox, part = 'gas'

    IF n_elements(cold) LT 5 THEN halfgr = 0 ELSE BEGIN
        histgradcum = weighted_histogram(sqrt(g[cold].x*g[cold].x + g[cold].y*g[cold].y)*lunit,nbins = 1000,max = stat[statind].rvir,weight = g[cold].mass,locations = locations,/cum)
        IF histgradcum[999]/total(g[cold].mass) LT 0.5 THEN STOP
        uniqind = uniq(histgradcum)
        locations = locations[uniqind]
        histgradcum = histgradcum[uniqind]
        temp = min(abs(histgradcum/total(g[cold].mass) - 0.5),ind5)
        IF (ind5 LT 2) THEN stop
        halfgr = spline(histgradcum[ind5 - 2:ind5 + 2]/total(g[cold].mass),locations[ind5 - 2:ind5 + 2],0.5)
;        halfgr = spline(histgradcum[uniqind]/total(g[cold].mass),locations[uniqind],0.5)
        IF NOT finite(halfgr) OR halfgr LT 0 THEN stop
    ENDELSE

    IF n_elements(s) LT 5 THEN halfsr = 0 ELSE BEGIN
        histsradcum = weighted_histogram(sqrt(s.x*s.x + s.y*s.y)*lunit,nbins = 1000,max = stat[statind].rvir,weight = s.mass,locations = locations,/cum)
        IF histsradcum[999]/total(s.mass) LT 0.5 THEN STOP
        uniqind = uniq(histsradcum)
        locations = locations[uniqind]
        histsradcum = histsradcum[uniqind]
        temp = min(abs(histsradcum/total(s.mass) - 0.5),ind5)
        IF (ind5 LT 2) THEN stop
        halfsr = spline(histsradcum[ind5 - 2:ind5 + 2]/total(s.mass),locations[ind5 - 2:ind5 + 2],0.5)
        IF NOT finite(halfsr) OR halfsr LT 0 THEN stop
   ENDELSE

    openw,1,fileroot + '.' + haloid + '_prof.txt',/APPEND
    printf,1,filebase,validids[i],1.0/h.time - 1.0,stat[statind].vc,total([d.mass,s.mass,g.mass])*munit,total(d.mass)*munit,total(s.mass)*munit,total(g.mass)*munit,gmass_HI,gmass_H2,gmass_cold,gmass_hot,gmassz,gmass_coldz,gmass_hotz,halfgr,halfsr,format='(A40,16A28)'
 ;   printf,1,'Gas Mass [Msun]','HI Mass [Msun]','H2 Mass [Msun]','Half Cold Gas Radius [kpc]','Half Cold Gas Height [kpc]','Log Gas Density Distro Mean ','Log Gas Density Distro STDEV ','Total Gas Metallicity ','Ox ','Fe ','SFR ',format='(11A40)'
;    printf,1,total(g.mass),total(g.mass*HI),total(g.mass*H2),g_halflengthall,g_halfheightall,pdfterms[1],pdfterms[2],total(g.mass*sfe*(2.09*Ox + 1.06*Fe))/total(g.mass*sfe),total(g.mass*sfe*ox)/total(g.mass*sfe),total(g.mass*sfe*fe)/total(g.mass*sfe),total(smassform[currentSF])/dtime,format='(11A40)'
    close,1
    print,filebase,validids[i],1.0/h.time - 1.0,stat[statind].vc,total([d.mass,s.mass,g.mass])*munit,total(d.mass)*munit,total(s.mass)*munit,total(g.mass)*munit,gmass_HI,gmass_H2,gmass_cold,gmass_hot,gmassz,gmass_coldz,gmass_hotz,halfgr,halfsr,format='(A40,16A18)'
;    stop
ENDFOR
close,/all
;stop
END
