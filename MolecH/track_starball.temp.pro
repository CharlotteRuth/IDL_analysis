PRO track_starball_master,outplot = outplot
!Y.STYLE = 1
!X.STYLE = 1
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    !X.MARGIN = [10,3]
    !Y.MARGIN = [4,2]
ENDELSE

pfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/h516.cosmo25cmb.2304g14HBWK.param'
units = tipsyunits(pfile)
msol_per_sysmass       = units.MASSUNIT ;2.310e15
kpc_per_syslength        = units.lengthunit ;25000.

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

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
;file = "h516.cosmo25cmb.1536g14HBWK"
;steps = ['00024','00048','00072','00092','00096','00100','00103','00104','00108','00112','00116','00120','00128','00132','00136','00140','00144','00148','00152','00156','00160','00164','00167','00168','00172','00176','00180','00184','00188','00192','00196','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
;file = "h516.cosmo25cmb.1536g3HBWK"
;steps = ['00030','00045','00060','00075','00090','00105','00120','00135','00150','00165','00180','00195','00210','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00468','00480','00492','00504','00512']

;file = "h516.cosmo25cmb.1536g1MBWK"
;steps = ['00024','00037','00046','00048','00061','00072','00084','00096','00120','00128','00144','00168','00192','00216','00227','00240','00264','00271','00288','00312','00328','00336','00360','00384','00406','00408','00432','00455','00456','00480','00504','00512']

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.2304g/'
file = 'h516.cosmo25cmb.2304g14HBWK'
;steps =
;['00024','00036','00037','00046','00048','00060','00061','00072','00084','00096','00108','00120','00132','00144','00148','00160','00172','00184',
;steps
;=['00196','00212','00224','00232','00240','00244','00256','00268','00272','00276','00288','00300','00308','00312','00324','00336','00348','00360','00372','00384','00396','00406','00416','00428','00440','00448']
;steps =['00196','00200','00204','00208','00212','00216','00220','00224','00227','00228','00232','00236','00240','00244','00248','00252','00256','00260','00264','00268','00271','00272','00276','00280','00284','00288','00292','00296','00300','00304','00308','00312','00316','00320','00324','00328','00332','00336','00340','00344','00348','00352','00356','00360','00364','00368','00372','00376','00380','00384','00388','00392','00396','00400','00404','00406','00408','00412','00416','00420','00424','00428','00432','00436','00440','00444','00448','00452','00455','00456','00460','00464','00468','00472',
;steps = ['00276','00388','00476','00480','00484','00488','00492','00496','00500','00504','00508','00512']
steps = ['00248','00252','00256','00260','00264','00268','00271','00272','00276','00512']


dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/'
file = 'h516.cosmo25cmb.3072g1MBWK'
steps = ['00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00492']
halo = ['1','1','1','1','1','1','1','1','1','1','1']
closerange = 5 
farrange = 10

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/'
file = 'h799.cosmo25cmb.3072g1MBWK'
;steps =['00024','00036','00048','00060','00072','00084','00096','00108',
steps =['00120','00132','00144','00156','00168','00180','00192','00204','00512']
;halo = ['1','1','1','1','1','2','1','2',
halo = ['2','1','1','1','1','1','1']
closerange = 2 
farrange = 5

sfile = dir + file + '/' + file + '.starlog'
mass = fltarr(N_ELEMENTS(steps),4)
mfile = dir + file + '/steps/'+ file + '.' +steps[N_ELEMENTS(steps) - 1] + '.dir/' + file + '.' +steps[N_ELEMENTS(steps) - 1] ;+ '.halo.1.starball.mrk'
FOR i = 0,N_ELEMENTS(steps) - 1 DO BEGIN
;    dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.'+step+'.dir/'
;    tfile = 'h516.cosmo25cmb.1536g3HBWK.'+step
;    sfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK.starlog'
;    mdir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/'

;    dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noJeans/h516.cosmo25cmb.1536g3HBWK_noJeans.'+step+'.dir/'
;    tfile = 'h516.cosmo25cmb.1536g3HBWK_noJeans.'+step
;    sfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK_noJeans.starlog'
;    mdir = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noJeans/h516.cosmo25cmb.1536g3HBWK_noJeans.00512.dir/' 

    tfile = dir + file + '/steps/'+ file + '.' +steps[i] + '.dir/' + file + '.' +steps[i]
    track_starball, tfile, sfile, msol_per_sysmass, kpc_per_syslength, name = file + '.' +steps[i], /molecularH, mfile = mfile,outplot = outplot,halo = halo[i],farrange = farrange,closerange = closerange;,/molecularH 

    data = read_stat_struc_amiga(dir + file + '/steps/'+ file + '.' +steps[i] + '.dir/' + file + '.' +steps[i] + '.amiga.stat')
    mass[i,0] = data[0].m_tot
    mass[i,1] = data[0].m_gas
    mass[i,2] = data[0].m_star
    mass[i,3] = data[0].m_tot - data[0].m_gas - data[0].m_star
    ind = where(mass[*,0] gt 0)
    IF N_ELEMENTS(ind) gt 1 THEN BEGIN
        IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename=tfile+'_mass.eps',/times,ysize=5,xsize=8,/inch ELSE window,3,xsize = 630,ysize = 392  
        !p.multi = [0,2,1]
        plot,steps, mass[ind,3], psym = -2, xstyle = 1, xrange = [0,512], xtitle = 'Step', ytitle = 'Mass [Msol]',yrange = [0,4e10]
        oplot,steps,mass[ind,1] + mass[*,2],psym = -1
        legend,['Dark Matter','Baryonic Mass'],psym = [-2,-1]
        
        plot,steps,mass[ind,1],psym = -4,xstyle = 1, xrange = [0,512], xtitle = 'Step', ytitle = 'Mass [Msol]',yrange = [0,4e9]
        oplot,steps,mass[ind,2],psym = -5
        legend,['Gas','Stars'],psym = [-4,-5],/right
        !p.multi = 0 

        IF KEYWORD_SET(outplot) THEN device,/close
    ENDIF
ENDFOR
END

PRO track_starball, tfile, sfile, msol_per_sysmass, kpc_per_syslength_a, max_tform = max_tform, mfile = mfile, name = name, outplot = outplot, MOLECULARH = MOLECULARH,halo = halo,closerange = closerange,farrange = farrange;mfile_halo = mfile_halo,
!Y.STYLE = 1
!X.STYLE = 1
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.0;2.25
    !y.charsize=1.0;2.25
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.0
    !y.charsize=1.0  
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
ENDELSE

hubble = 73.0
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
dDelta = 0.003
sec_per_year = 31556926
molec_weight = (0.76*1 + 0.24*4.0)
densunit = msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength_a^3/cm_per_kpc^3
timeunit=SQRT((kpc_per_syslength_a*3.086d21)^3/(6.67d-8*msol_per_sysmass*1.99d33))/(3600.*24.*365.24)

IF NOT (keyword_set(halo)) THEN halo = '1'
IF NOT (keyword_set(closerange)) THEN closerange = 5
IF NOT (keyword_set(farrange)) THEN farrange = 25
 
print,tfile
rtipsy,tfile,header,g,d,s,/JUSTHEAD
rtipsy,mfile,header_mrk,g,d,s,/JUSTHEAD
rtipsy,tfile+'.halo.'+halo,h,g,d,s 
;ntot = header.n
ndark = header.ndark
;nstar = h.nstar
;ngas = h.ngas
print,'ntot: ',strtrim(h.n,2),', ndark ',strtrim(h.ndark,2),', ngas: ',strtrim(h.ngas,2),', nstar: ',strtrim(h.nstar,2)
IF KEYWORD_SET(mfile) THEN BEGIN
    readarr,tfile+'.halo.'+halo+'.iord',h,iord_gas,type = 'int',/ascii,part = 'gas'
    readarr,tfile+'.halo.'+halo+'.iord',h,iord_star,type = 'int',/ascii,part = 'star'
    readarr,mfile+'.iord',header_mrk,iord_all_all,type = 'long',/ascii
ENDIF

kpc_per_syslength = kpc_per_syslength_a*h.time
s.x = s.x*kpc_per_syslength
s.y = s.y*kpc_per_syslength
s.z = s.z*kpc_per_syslength
d.x = d.x*kpc_per_syslength
d.y = d.y*kpc_per_syslength
d.z = d.z*kpc_per_syslength
g.x = g.x*kpc_per_syslength
g.y = g.y*kpc_per_syslength
g.z = g.z*kpc_per_syslength

star_info = rstarlog(sfile);,MOLECULARH = MOLECULARH)
IF NOT KEYWORD_SET(max_tform) THEN max_tform = MAX(star_info.timeform)
star_info = star_info[where(star_info.iorderstar ne 0)]

IF KEYWORD_SET(mfile) THEN BEGIN
    openr,1,mfile+ '.halo.1.starball.mrk'
    readf,1,ntot,ngas,nstar
    close,1
    readcol,mfile+ '.halo.1.starball.mrk',mark
    mrkind = mark[1:N_ELEMENTS(mark)-1] - 1

    siord_all = iord_all_all[mrkind[where(mrkind gt ndark + ngas AND mrkind le ndark + ngas + nstar)]]
    sl_starball_ind = [0]
    FOR is = 0, N_ELEMENTS(siord_all) - 1 DO BEGIN
        ind = (where(star_info.iorderstar eq siord_all[is]))[0]
        sl_starball_ind = [sl_starball_ind,ind]
    ENDFOR
    sl_starball = star_info[sl_starball_ind[1:N_ELEMENTS(sl_starball_ind) - 1]]
    IF (MIN(siord_all) LT MAX(iord_star)) THEN BEGIN
        siord= siord_all[where(siord_all le max(iord_star))]
        sind = fltarr(N_ELEMENTS(siord))
        FOR is = 0, N_ELEMENTS(sind) - 1 DO BEGIN
            sind[is] = where(siord_all[is] eq iord_star)
        ENDFOR
        s_starball = s[sind]
        nostar = 0
    ENDIF ELSE BEGIN
        sind = [0]
        nostar = 1
    ENDELSE

    giord = iord_all_all[mrkind[where(mrkind le ngas)]]
    IF (giord[0] ne -1) THEN BEGIN
        gind = fltarr(N_ELEMENTS(giord))
        FOR ig = 0, N_ELEMENTS(giord) - 1 DO BEGIN
            IF (where(giord[ig] eq iord_gas))[0] NE -1 THEN gind[ig] = (where(giord[ig] eq iord_gas))[0]
        ENDFOR
    ENDIF
    g_starball= g[gind]
    gprior_starball = [0]
    g_timeform = [0]
    print,'# of star in ball: ',STRTRIM(N_ELEMENTS(siord_all),2),', # of stars formed so far: ',STRTRIM(N_ELEMENTS(siord),2)
    FOR is = N_ELEMENTS(siord), N_ELEMENTS(siord_all) - 1 DO BEGIN
        ind = (where(star_info.iorderstar eq siord_all[is]))[0]
        if ind[0] ne -1 THEN BEGIN
            sfiord = star_info[ind].IORDERGAS
            gind = where(sfiord eq iord_gas)
            gprior_starball = [gprior_starball,gind]
            g_timeform = [g_timeform,star_info[ind].timeform]
        ENDIF ELSE print,'no ind'
        if is MOD 500 eq 0 then print,is
    ENDFOR
    IF (N_ELEMENTS(gprior_starball) gt 1) THEN gprior_starball = gprior_starball[1:N_ELEMENTS(gprior_starball) - 1]
    g_starball = [g[gprior_starball],g_starball]
    blank = fltarr(N_ELEMENTS(g_starball))
    IF (N_ELEMENTS(gprior_starball) gt 1) THEN g_timeform = [g_timeform[1:N_ELEMENTS(g_timeform) - 1],blank] ELSE g_timeform = blank
ENDIF

IF KEYWORD_SET(mfile_halo) THEN BEGIN
    openr,1,mfile_halo
    readf,1,ntot,ngas,nstar
    close,1
    readcol,mfile_halo,mark
    iord = mark[1:N_ELEMENTS(mark)-1] - 1
    giord = iord[where(iord le ngas)]
    diord = iord[where(iord gt ngas AND iord le ndark + ngas)] - ngas
    siord_all = iord[where(iord gt ndark + ngas AND iord le ndark + ngas + nstar)] - ndark - ngas
    siord = siord_all[where(siord_all le h.nstar)]
    g_halo= g[giord]
    d_halo = d[diord]
    s_halo = s[siord]
;    IF N_ELEMENTS(s_halo) gt 500 THEN BEGIN
;        x_center = MEAN(s_halo[0:500].x)
;        y_center = MEAN(s_halo[0:500].y)
;    ENDIF ELSE BEGIN

ENDIF ELSE BEGIN
    g_halo= g
    d_halo = d
    s_halo = s
    x_center = 0
    y_center = 0     
    z_center = 0   
ENDELSE
s_halo.x = s_halo.x - x_center
s_halo.y = s_halo.y - y_center
s_halo.z = s_halo.z - z_center
d_halo.x = d_halo.x - x_center
d_halo.y = d_halo.y - y_center
d_halo.z = d_halo.z - z_center
g_halo.x = g_halo.x - x_center
g_halo.y = g_halo.y - y_center
g_halo.z = g_halo.z - z_center


;---------------------------------------- Frame -----------------------------------
;!p.multi = [0,2,2]
!X.MARGIN = [6,1]
!Y.MARGIN = [3,2]
IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename = tfile+'.halo.1_image.eps',/times,ysize=18,xsize=18 ELSE window,0,xsize = 800,ysize = 800 ;392
multiplot,[2,2],gap=0.04,mtitle = name,/square,/doxaxis,/doyaxis
dxy = closerange           ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
loadct,0
plot,s_halo.x,s_halo.y,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],ytitle = 'Y [kpc]';,title = name
oplot,g_halo.x,g_halo.y,psym = 3,color = 150
oplot,s_halo.x,s_halo.y,psym = 3
IF KEYWORD_SET(mfile) THEN BEGIN
    loadct,18
;oplot,g_starball.x,g_starball.y,psym = 3,color = 244
;g_color = FIX(g_timeform/max_tform*220.0) + 24
;g_color =  FIX((g_timeform - h.time)*170.0)
    ind = g_timeform[where(g_timeform eq 0)]
    g_color =  FIX((g_timeform - MIN(g_timeform))/max_tform*129.0)
    g_color[ind] = 129
    FOR ig = 0,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].y,g_starball[ig].y],psym = 3,color = g_color[ig]
    IF NOT nostar THEN BEGIN
;s_color = FIX(s_starball.tform/max_tform*220.0) + 24
;s_color = 299 - FIX(s_starball.tform/max_tform*170.0)
        s_color = 130 + FIX((MAX(s_starball.tform) - s_starball.tform)/max_tform*130.0)
        FOR is = 0,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].y,s_starball[is].y],psym = 3,color = s_color[is]
    ENDIF
ENDIF ELSE BEGIN
    loadct,18
    s_color = 130 + FIX((MAX(s_halo.tform) - s_halo.tform)/max_tform*299.0)
    FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].y,s_halo[is].y],psym = 3,color = s_color[is]
ENDELSE       
multiplot,/doxaxis,/doyaxis

dxy = farrange          ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
loadct,0
plot,s_halo.x,s_halo.y,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy];,title = name
oplot,g_halo.x,g_halo.y,psym = 3,color = 150
oplot,s_halo.x,s_halo.y,psym = 3
IF KEYWORD_SET(mfile) THEN BEGIN
    loadct,18
    FOR ig = 0,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].y,g_starball[ig].y],psym = 3,color = g_color[ig]
    IF NOT nostar THEN BEGIN
        FOR is = 0,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].y,s_starball[is].y],psym = 3,color = s_color[is]
    ENDIF
ENDIF ELSE BEGIN
    loadct,18
    FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].y,s_halo[is].y],psym = 3,color = s_color[is]
ENDELSE 
multiplot,/doxaxis,/doyaxis    

dxy = closerange           ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
loadct,0
plot,s_halo.x,s_halo.z,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],xtitle = 'X [kpc]',ytitle = 'Z [kpc]';,title = name
oplot,g_halo.x,g_halo.z,psym = 3,color = 150
oplot,s_halo.x,s_halo.z,psym = 3
IF KEYWORD_SET(mfile) THEN BEGIN
    loadct,18
    FOR ig = 0,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].z,g_starball[ig].z],psym = 3,color = g_color[ig]
    IF NOT nostar THEN BEGIN
        FOR is = 0,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].z,s_starball[is].z],psym = 3,color = s_color[is]
    ENDIF
ENDIF ELSE BEGIN
    loadct,18
    FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].z,s_halo[is].z],psym = 3,color = s_color[is]
ENDELSE 
multiplot,/doxaxis,/doyaxis
    
dxy = farrange          ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
loadct,0
plot,s_halo.x,s_halo.z,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],xtitle = 'X [kpc]';,title = name
oplot,g_halo.x,g_halo.z,psym = 3,color = 150
oplot,s_halo.x,s_halo.z,psym = 3
IF KEYWORD_SET(mfile) THEN BEGIN
    loadct,18
    FOR ig = 0,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].z,g_starball[ig].z],psym = 3,color = g_color[ig]
    IF NOT nostar THEN BEGIN
        FOR is = 0,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].z,s_starball[is].z],psym = 3,color = s_color[is]
    ENDIF
ENDIF ELSE BEGIN
    loadct,18
    FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].z,s_halo[is].z],psym = 3,color = s_color[is]
ENDELSE 

IF KEYWORD_SET(outplot) THEN  device,/close
multiplot,/reset

;---------------------- SFH ------------------------------------------------------------------------
loadct,39
!p.multi = 0
!X.MARGIN = [12,3]
!Y.MARGIN = [6,2]
IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename=tfile+'.halo.1_sfh.eps',/times,ysize=12,xsize=18 ELSE window,2,xsize = 630,ysize = 392
sfr,s_halo,massunit = msol_per_sysmass,timeunit = timeunit,binsize = 1e7,xrange = [0,5.5],yrange = [0,0.25],xrnage = [0,timeunit*max_tform/1e9],massform = MAX(s.mass)
IF (KEYWORD_SET(mfile)) THEN $
  IF NOT nostar THEN sfr,s_starball,massunit = msol_per_sysmass,timeunit = timeunit,/overplot,color = 240,binsize = 1e7
IF KEYWORD_SET(outplot) THEN device,/close else stop
END
