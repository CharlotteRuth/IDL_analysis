PRO track_mass_master

spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

;dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
;file = 'h516.cosmo25cmb.1536g6MbwK'
;steps = ['00030','00045','00060','00075',$
;         '00090','00105','00120','00135','00150',$
;         '00165','00180','00195','00210','00216','00228','00240','00252',$
;         '00264','00276','00288','00300','00312','00324','00336',$
;       '00348','00360','00372','00384','00396','00408','00420',$
;         '00432','00444','00456','00468','00480','00492','00504','00512']
;['00105','00120','00135','00150','00195','00252','00300','00348','00396','00456','00504','00512']
;['00075','00090','00105','00120','00135','00150','00195','00252','00300','00348','00396','00456','00504','00512']
;steps = ['00024','00048','00060','00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00456','00468','00480','00492','00504','00512']
;steps = ['00192','00240','00328','00384','00406','00455','00512']

;dir = prefix + 'h516.cosmo25cmb.1536g/'
;file = "h516.cosmo25cmb.1536g14HBWK"
;steps = ['00024','00048','00072','00092','00096','00100','00103','00104','00108','00112','00116','00120','00128','00132','00136','00140','00144','00148','00152','00156','00160','00164','00167','00168','00172','00176','00180','00184','00188','00192','00196','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']

;dir = prefix + 'h516.cosmo25cmb.1536g/'
;file = "h516.cosmo25cmb.1536g3HBWK"
;steps = ['00030','00045','00060','00075','00090','00105','00120','00135','00150','00165','00180','00195','00210','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00468','00480','00492','00504','00512']

;file = "h516.cosmo25cmb.1536g1MBWK"
;steps = ['00024','00037','00046','00048','00061','00072','00084','00096','00120','00128','00144','00168','00192','00216','00227','00240','00264','00271','00288','00312','00328','00336','00360','00384','00406','00408','00432','00455','00456','00480','00504','00512']

;dir = prefix + 'h516.cosmo25cmb.2304g/'
;file = 'h516.cosmo25cmb.2304g14HBWK'
;steps =
;['00024','00036','00037','00046','00048','00060','00061','00072','00084','00096','00108','00120','00132','00144','00148','00160','00172','00184',
;steps
;=['00196','00212','00224','00232','00240','00244','00256','00268','00272','00276','00288','00300','00308','00312','00324','00336','00348','00360','00372','00384','00396','00406','00416','00428','00440','00448']
;steps =['00196','00200','00204','00208','00212','00216','00220','00224','00227','00228','00232','00236','00240','00244','00248','00252','00256','00260','00264','00268','00271','00272','00276','00280','00284','00288','00292','00296','00300','00304','00308','00312','00316','00320','00324','00328','00332','00336','00340','00344','00348','00352','00356','00360','00364','00368','00372','00376','00380','00384','00388','00392','00396','00400','00404','00406','00408','00412','00416','00420','00424','00428','00432','00436','00440','00444','00448','00452','00455','00456','00460','00464','00468','00472',
;steps = ['00276','00388','00476','00480','00484','00488','00492','00496','00500','00504','00508','00512']
;steps = ['00248','00252','00256','00260','00264','00268','00271','00272','00276','00512']

;dir = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/'
;file = 'h516.cosmo25cmb.3072g1MBWK'
;steps_st = ['00012','00024','00036','00048','00060','00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00456','00468','00480','00492','00504','00512']
;halo = ['2','2','2','3','2','2','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

;dir = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/'
;file = 'h516.cosmo25cmb.3072g14HBWK'
;steps_st = ['00024','00032','00080','00088','00092','00104','00128','00132','00144','00148','00156','00168','00180','00192','00204','00216','00228','00240','00252','00260','00272','00284','00308','00320','00332','00344','00356','00368','00380','00392','00404','00420','00432','00436','00444','00456','00468','00480','00492','00504','00512']
;halo = ['2','4','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']
;steps_st = ['00512']
;halo = ['1']

;dir = prefix + 'h799.cosmo25cmb.3072g/'
;file = 'h799.cosmo25cmb.3072g1MBWK'
;steps =['00024','00036','00048','00060','00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00512']
;halo = ['1',    '1',    '1',    '1',    '1',    '2',    '1',    '2',    '2',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1']

;dir = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/'
;file = 'h799.cosmo25cmb.3072g14HBWK'
;steps_st =['00048','00060','00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00456','00468','00480','00492','00504','00512']
;halo =    ['1',    '1',    '1',    '2',    '2',    '2',    '2',    '4',    '4',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1']

;dir = prefix + 'h603.cosmo50cmb.3072g/'
;file = 'h603.cosmo50cmb.3072gs1MbwK'
;steps = ['00348','00420','00432','00444','00468','00480','00492','00504','00512']
;halo = ['1','1','1','1','1','1','1','1','1']

;dir = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/'
;file = 'h603.cosmo50cmb.3072g14HBWK'
;steps_st = ['00072','00084','00088','00092','00096','00100','00104','00108','00120','00132','00144','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00456','00468','00480','00492','00504','00512']
;halo =  ['2','2','2','2','2','2','2','2','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

;dir = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/'
;file = 'h603.cosmo50cmb.3072gs1MbwK'
;steps_st = ['00036','00048','00060','00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00456','00468','00480','00492','00504','00512']
;halo = ['2','2','3','2','7','2','2','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']

dir = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/'
file = 'h986.cosmo50cmb.3072g14HBWK'
steps_st = ['00012','00024','00036','00048','00060','00072','00084','00092','00108','00120','00132','00144','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00396','00408','00432','00444','00456','00468','00492','00504','00512']
halo     = [    '1',    '1',    '1',    '1',    '2',    '2',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1']

j = track_j(dir,file,steps_st = steps_st,halo = halo)
writecol,'j.dat',j
stop
mass = track_mass(dir,file,steps_st = steps_st,halo = halo,coldgas = 1)
openw,1,dir + 'masstrack.dat'
printf,1,transpose(mass),format = '('+strtrim((size(mass))[2],2)+'E)'
close,1
stop
track_mass_metal,dir,file,steps_st = steps_st,halo = halo
END

PRO track_mass_plot,dir,file,outplot = outplot,key = keys,linestyle = linestyles,thick = thicks,colors = colors,steps_st = steps_st

!Y.STYLE = 1
!X.STYLE = 1
n = N_ELEMENTS(file)
formatplot,outplot = outplot
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    if colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    IF NOT KEYWORD_SET(thick) THEN thick = findgen(n)*2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = REVERSE(findgen(n)*2   )
ENDELSE

FOR i = 0, n - 1 DO BEGIN
    dir1 = dir[i]
    cut_pos = strsplit(dir1,'/')
    dir1_cut = strmid(dir1,0,cut_pos[N_ELEMENTS(cut_pos) - 1])

    file1 = file[i]
    cut_pos = strsplit(file1,'.')
    file1_cut = strmid(file1,0,cut_pos[N_ELEMENTS(cut_pos) - 1] - 1)
    print,dir1_cut,file1_cut
    mass = track_mass(dir1_cut,file1_cut,coldgas = 1)
    
    openw,1,dir1_cut+'masstrack.dat'
    printf,1,transpose(mass),format = '('+strtrim((size(mass))[2],2)+'E)'
    close,1
ENDFOR

END


FUNCTION track_mass,dir,file,steps_st = steps_st,halo = halo,outplot = outplot,coldgas = coldgas
formatplot,outplot = outplot

grav = 4.51737015d-39;kpc^3 msol^(-1) s^(-2)
kB = 1.3806503d-16 ;cm^2 g s^-2 K^-1
g_per_msol = 1.98892d33
cm_per_kpc = 3.08568025d21 

IF NOT KEYWORD_SET(steps_st) THEN BEGIN
;   command = 'ls ' + dir + file + "/steps/"+"*.dir/*amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
    command = 'ls ' + dir + "*.dir/*amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
;   command = 'ls ' + dir + "*/*halo.1.std | cut -d '.' -f7"
    spawn,command,files       ; return the files without the .amiga.grp
    grpfile = files+'.amiga.grp'
    statfile = files+'.amiga.stat'
    file_temp = files[0]
    cut_pos = strsplit(file_temp,'.')
    files_cut = strmid(file_temp,0)
    steps_st = strmid(files,cut_pos[N_ELEMENTS(cut_pos) - 1])
;    step_st = strmid(files,4,5,/reverse_offset)
;    step = fix(strmid(files,2,3,/reverse_offset))
;    stop
ENDIF 
steps = fix(steps_st)

IF NOT KEYWORD_SET(halo) THEN BEGIN
    halo = strarr(N_ELEMENTS(steps)) + '1'
ENDIF
IF KEYWORD_SET(coldgas) THEN mass = fltarr(N_ELEMENTS(steps),14)  ELSE mass = fltarr(N_ELEMENTS(steps),7)

;readcol,dir + file + '/' + file + '.log',dTime, z, E, T, U, Eth, Lx, Ly, Lz, WallTime, dWMax, dImax, dEMax, dMultiEff,/silent
readcol,dir + '../' + file + '.log',dTime, z, E, T, U, Eth, Lx, Ly, Lz, WallTime, dWMax, dImax, dEMax, dMultiEff,/silent
units = tipsyunits(dir + '../' + file + '.param')
time0 = 0
FOR i = 0,N_ELEMENTS(steps_st) - 1 DO BEGIN
;    tfile = dir + file + '/steps_st/'+ file + '.' +steps_st[i] +
;    '.dir/' + file + '.' +steps_st[i]
    tfile = dir +  file + '.' +steps_st[i] + '.dir/' + file + '.' +steps_st[i]; + '.halo.'+halo
    IF KEYWORD_SET(coldgas) THEN BEGIN
        rtipsy,tfile + '.halo.'+halo[i] + '.std',h,g,d,s;,/justhead
        temp = min(abs(z - (1.0/h.time - 1.0)),ind)
        time = dTime[ind]*units.timeunit
        readarr,tfile + '.halo.'+halo[i]+".HI",h,HI,/ascii,part = 'gas'
;        readarr,tfile +  '.halo.'+halo[i]+".amiga.grp",h,grp,/ascii,part = 'gas' removed for h986
        IF (FILE_TEST(tfile+ '.halo.'+halo[i]+".H2")) THEN  readarr,tfile+ '.halo.'+halo[i]+".H2",h,H2,/ascii,part = 'gas' ELSE H2 = HI * 0
        H2 = H2*2
        coldg = TOTAL(H2*g.mass + HI*g.mass)
        coldg = coldg*units.massunit
;        coldg = TOTAL(H2 + HI)*units.gasparmass
        metal = TOTAL(H2*g.mass*g.zmetal + HI*g.mass*g.zmetal)/TOTAL(H2*g.mass + HI*g.mass)
;        metal = 0
        H2mass = TOTAL(H2*g.mass)*units.massunit
;        H2mass = TOTAL(H2)*units.gasparmass
        if h.nstar ne 0 THEN SFR = N_ELEMENTS(where(s.tform*units.timeunit gt time0))*units.istarmass/(time - time0) ELSE SFR = 0
;----------------------- Radius  -- The radii that contians half the
;                        disk gas particles
        g.x = g.x*units.lengthunit*h.time
        g.y = g.y*units.lengthunit*h.time
        g.z = g.z*units.lengthunit*h.time
        if h.nstar ne 0 THEN BEGIN
            s.x = s.x*units.lengthunit*h.time
            s.y = s.y*units.lengthunit*h.time
            s.z = s.z*units.lengthunit*h.time
        ENDIF
        radius_g = sqrt(g.x*g.x + g.y*g.y)
        ind_gall = where(ABS(g.z) lt 3 AND g.tempg le 1e4); AND grp eq FIX(halo[i])) removed for h986
        gall = g[ind_gall]
        radius_gall = radius_g[ind_gall]
        radius_gall = radius_gall[SORT(radius_gall)]
        r25 = radius_gall[N_ELEMENTS(radius_gall)/4]

;----------------------- Pressure
        gind = where(sqrt(g.x*g.x + g.y*g.y) le r25)
        sigma_g_const = 11 ;km/s estimate used in Bigiel 07
        sd_HI = TOTAL(HI[gind]*g[gind].mass)*units.massunit/(!PI*r25*r25)
        sd_H2 = TOTAL(H2[gind]*g[gind].mass)*units.massunit/(!PI*r25*r25)
        if h.nstar ne 0 THEN BEGIN 
            sind = where(sqrt(s.x*s.x + s.y*s.y) le r25)
            sv = sqrt(s[sind].vx*s[sind].vx + s[sind].vy*s[sind].vy + s[sind].vz*s[sind].vz)*units.vunit*h.time
            sigma_sx = sqrt((moment(s[sind].vx*units.vunit*h.time))[1])
            sigma_sy = sqrt((moment(s[sind].vy*units.vunit*h.time))[1])
            sigma_sz = sqrt((moment(s[sind].vz*units.vunit*h.time))[1]) ;line of sight velocity dispersion
            sigma_s = sqrt(sigma_sx*sigma_sx + sigma_sy*sigma_sy + sigma_sz*sigma_sz)
            sd_s =  TOTAL(         s[sind].mass)*units.massunit/(!PI*r25*r25)

            gind = where(sqrt(g.x*g.x + g.y*g.y) le r25)
            sv = sqrt(g[gind].vx*g[gind].vx + g[gind].vy*g[gind].vy + g[gind].vz*g[gind].vz)*units.vunit*h.time
            sigma_gx = sqrt((moment(g[gind].vx*units.vunit*h.time))[1])
            sigma_gy = sqrt((moment(g[gind].vy*units.vunit*h.time))[1])
            sigma_gz = sqrt((moment(g[gind].vz*units.vunit*h.time))[1]) ;line of sight velocity dispersion                 
            sigma_g = sqrt(sigma_gx*sigma_gx + sigma_gy*sigma_gy + sigma_gz*sigma_gz)

            pressure = !PI/2*Grav*(sd_HI + sd_H2)*((sd_HI + sd_H2) + sd_s*sigma_g_const/sigma_sz)
            print,h.time,sigma_sz,pressure/cm_per_kpc*g_per_msol/kB,!PI/2*Grav*(sd_HI + sd_H2)*((sd_HI + sd_H2) + sd_s*sigma_g_const/sigma_sz*h.time)/cm_per_kpc*g_per_msol/kB
;            stop
        ENDIF ELSE pressure = !PI/2*Grav*(sd_HI + sd_H2)*((sd_HI + sd_H2))
        pressure = pressure/cm_per_kpc*g_per_msol/kB
        H2frac_center = sd_H2/sd_HI

;--------------------------       
        mass[i,7] = metal
        mass[i,8] = coldg 
        mass[i,9] = H2mass
        mass[i,10] = H2frac_center
        mass[i,11] = pressure
        mass[i,12] = r25
        mass[i,13] = SFR
    ENDIF ELSE rtipsy,tfile,h,g,d,s,/justhead
    temp = min(abs(z - (1.0/h.time - 1.0)),ind)
    time = dTime[ind]*units.timeunit
    time0 = time
;    data = read_stat_struc_amiga(dir + file + '/steps/'+ file + '.' +steps_st[i] + '.dir/' + file + '.' +steps_st[i] + '.amiga.stat')
    data = read_stat_struc_amiga(dir + file + '.' +steps_st[i] + '.dir/' + file + '.' +steps_st[i] + '.amiga.stat')
    grp = where(data.group eq halo[i])
    mass[i,0] = FLOAT(halo[i])
    mass[i,1] = time/1e9
    mass[i,2] = z[i]
    mass[i,3] = data[grp].m_tot
    mass[i,4] = data[grp].m_gas
    mass[i,5] = data[grp].m_star
    mass[i,6] = data[grp].m_tot - data[0].m_gas - data[0].m_star
    print,steps_st[i],ind,h.time,(1.0/h.time - 1.0),z[ind],dTime[ind],data[grp].m_tot,data[grp].m_gas,data[grp].m_star,data[grp].m_tot - data[0].m_gas - data[0].m_star
    ind = where(mass[*,1] gt 0)
    IF N_ELEMENTS(ind) gt 1 THEN BEGIN
        IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename=tfile+'.halo.1_mass' + outplot +'.eps',/times,ysize=5,xsize=8,/inch ELSE window,3,xsize = 630,ysize = 392  
        !p.multi = 0
;       !p.multi = [0,2,1]
        plot,mass[ind,1], mass[ind,3], psym = -2, xstyle = 1, xrange = [0,14], xtitle = 'Time [Gyr]', ytitle = 'Mass [Msol]',yrange = [0,2.5e11] ;[0,2.5e10]
        oplot,mass[ind,1],mass[ind,3] + mass[ind,3],psym = -1
        legend,['Dark Matter','Baryonic Mass'],psym = [-2,-1]
        
;        plot,mass[ind,0],mass[ind,2],psym = -4,xstyle = 1, xrange = [0,14], xtitle = 'Time [Gyr]', ytitle = 'Mass [Msol]',yrange = [0,2.5e10] ;[0,2.5e9]
;        oplot,mass[ind,0],mass[ind,3],psym = -5
;        legend,['Gas','Stars'],psym = [-4,-5],/right
;        !p.multi = 0 
        
        IF KEYWORD_SET(outplot) THEN device,/close
    ENDIF
ENDFOR
return,mass
END


PRO track_mass_metal,dir,file,steps_st = steps_st,halo = halo

IF NOT KEYWORD_SET(steps_st) THEN BEGIN
    command = 'ls ' + dir + "*.dir/*amiga.grp | grep amiga.grp | sed 's/.amiga.grp//g'"
    spawn,command,files       ; return the files without the .amiga.grp
    grpfile = files+'.amiga.grp'
    statfile = files+'.amiga.stat'
    file_temp = files[0]
    cut_pos = strsplit(file_temp,'.')
    files_cut = strmid(file_temp,0)
    steps_st = strmid(files,cut_pos[N_ELEMENTS(cut_pos) - 1])
ENDIF 
steps = fix(steps_st)

readcol,dir+'masstrack.dat',halo,time,z,mtot,mgas,mstar,mdark,metal,mcoldg,mH2,H2frac,Pressure,r25,SFR

mox = metal*0
FOR i = 0,N_ELEMENTS(steps_st) - 1 DO BEGIN
   mfile = dir +  file + '.' +steps_st[i] + '.dir/' + file + '.' +steps_st[i] + '.metals.fits'
   metaldat = mrdfits(mfile,1)
   ind = where(fix(metaldat.grp) EQ fix(halo[i]))
   IF (metaldat[ind].ox_sfr EQ 0 OR finite(metaldat[ind].ox_sfr) EQ 0) THEN $
      IF (metaldat[ind].ox_inner EQ 0 OR finite(metaldat[ind].ox_inner) EQ 0) THEN $
         IF (metaldat[ind].ox_inneratom EQ 0 OR finite(metaldat[ind].ox_inneratom) EQ 0) THEN mox[i] = metaldat[ind].ox_cold ELSE mox[i] = metaldat[ind].ox_inneratom $
      ELSE mox[i] = metaldat[ind].ox_inner $
      ELSE mox[i] = metaldat[ind].ox_sfr
ENDFOR
stop
mass = [[halo],[time],[z],[mtot],[mgas],[mstar],[mdark],[metal],[mox],[mcoldg],[mH2],[H2frac],[Pressure],[r25],[SFR]]
print,transpose(mass)
format = '(I,E,E,E,E,E,E,E,E,E,E,E,E,E,E)'
openw,1,dir+'masstrack.dat'
printf,1,transpose(mass),format = format
close,1
;writecol,'masstrack.dat',halo,time,z,mtot,mgas,mstar,mdark,metal,mox,mcoldg,mH2,H2frac,Pressure,r25,SFR,format = '(I,I,D,D,D,D,D,D,D,D,D,D,D,D,D)'
stop
END

FUNCTION track_j,dir,file,steps_st = steps_st, halo = halo
s =indgen(n_elements(dir))*0.1
xsi = 1.0-1.25*(1.0-(1.25-1.0)*alog(1.25/(1.25-1.0)))
ps = xsi*1.25*(1.25-1.0)/(xsi*s+1.25-1.0)^2.
loadct,39
steps = fix(steps_st)
units = tipsyunits(dir + '../' + file + '.param')

n = n_elements(steps)
j = fltarr(n)
IF NOT KEYWORD_SET(halo) THEN BEGIN
 hhalo = strarr(n) + '1'
ENDIF

FOR i = 0, n -1 DO BEGIN
   cd,dir
;   ytickname=[' ','0.2',' ','0.6',' ','1.0'] 
;   jdist, prefix=dir+file + '.' + steps_st[i] + '.dir/', file + '.' + steps_st[i], halo = halo[i], h, hs, hg, d, s, g, fdisk, dmass, smass, gmass, /obs 
;   jdist_bary = (gmass/max(gmass)+smass/max(smass))*(fdisk/.2121)
;   plot, indgen(n_elements(h))*0.05,dmass/max(dmass), xrange=[0,2.5], yrange=[0,1], xtickname=xtickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)',/nodata,title = label
;   oplot,indgen(n_elements(h))*0.05,dmass/max(dmass), linestyle = 2
;   if n_elements(hs) lt n_elements(hg) then b = n_elements(hg) else b = n_elements(hs)
;   oplot, indgen(b)*0.05,jdist_bary, linestyle=0
;   j[i] = total(jdist_bary*(indgen(b)*0.05))/total(jdist_bary)

   filename = dir+file+'.'+steps_st[i]+'.dir/'+file+'.'+steps_st[i]+'.halo.'+halo[i]
   rtipsy,filename,head,g,d,s
   readarr, filename + '.HI',head,hi,/ascii,part = 'gas'
   IF (file_test(filename + ".H2")) THEN  readarr,filename+".H2",head,H2,/ascii,part = 'gas' ELSE H2 = HI * 0
   angmom, g, jvecg, lvec, jmagg, ltot
   ind = where(jvecg[*,2] ge 0.)
   j[i] = total(jvecg[ind,2]*(2.0*H2 + HI)*head.time*head.time)/total(2.0*H2 + HI) ;distance*mass*speed
 ;  g = g[ind]
;   hg = histogram(jvecg[ind,2]/mean(jvecg[ind,2]), binsize=0.05, reverse_indices=ri)
    print,steps_st[i],j[i]
;   oplot,[j[i],j[i]],[0,1]
;   stop
ENDFOR
return,j
END

;h516.cosmo25cmb.3072g14HBWK.00132.dir/h516.cosmo25cmb.3072g14HBWK.00132.amiga.grp
;h516.cosmo25cmb.3072g14HBWK.00168.dir/h516.cosmo25cmb.3072g14HBWK.00168.amiga.grp
;h516.cosmo25cmb.3072g14HBWK.00344.dir/h516.cosmo25cmb.3072g14HBWK.00344.amiga.grp
;h516.cosmo25cmb.3072g14HBWK.00420.dir/h516.cosmo25cmb.3072g14HBWK.00420.amiga.grp
;h516.cosmo25cmb.3072g14HBWK.00468.dir/h516.cosmo25cmb.3072g14HBWK.00468.amiga.grp
;h516.cosmo25cmb.3072g14HBWK.00480.dir/h516.cosmo25cmb.3072g14HBWK.00480.amiga.grp
;h516.cosmo25cmb.3072g14HBWK.00492.dir/h516.cosmo25cmb.3072g14HBWK.00492.amiga.grp
;h516.cosmo25cmb.3072g14HBWK.00504.dir/h516.cosmo25cmb.3072g14HBWK.00504.amiga.grp
