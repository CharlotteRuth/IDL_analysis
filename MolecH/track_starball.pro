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

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/'
file = 'h603.cosmo50cmb.3072gs1MbwK'
steps = ['00348','00408','00420','00432','00444','00456','00468','00480','00492','00504','00512']
halo = ['1','1','1','1','1','1','1','1','1','1','1']
closerange = 2.5;10
farrange = 10;
molecularH = 0

;168
dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/'
file = 'h986.cosmo50cmb.3072g14HBWK'
steps =['00024','00036','00048','00060','00072','00084','00156','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00396','00408','00432','00444','00456','00468','00492','00504','00512']
halo = ['1','1','1','1','1','2','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']
closerange = 2.5;10
farrange = 10;
molecularH = 1

;decomp halo needed:'00132','00168' ,'00204' ,'00264','00288','00420'
;halo needed: '00096','00252','00348','00360'
;'2','1','1','1','1','1'

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/'
file = 'h603.cosmo50cmb.3072g14HBWK'
;steps =['00072','00084','00088','00092','00100','00104','00108','00120','00144','00156','00180','00192','00216','00228','00240','00264','00276','00300','00312','00324','00336','00372','00384','00396','00408','00432','00444','00456','00468','00480','00492','00504','00512']
;halo = ['2','2','2','2','2','2','2','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']
steps = ['00408','00432','00444','00456','00468','00480','00492','00504','00512']
halo = ['1','1','1','1','1','1','1','1','1']
;closerange = 2.5;10
;farrange = 10;
;molecularH = 1

colorarr_ext = 'decomp'
compext = 'thdisk';'thdisk';'bulge','disk','pb','halo'
mfext = '.decomp.'+ compext +'.mrk' ;create with select_comp

dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/'
file = 'h516.cosmo25cmb.3072g14HBWK'
steps = ['00228','00232','00236','00240','00244','00248','00252','00260','00264','00268','00272','00276','00280','00284','00308','00320','00332','00512']
;,'00344','00356','00368','00380','00392','00404','00420','00432','00436','00444','00456','00468','00480','00492','00504','00512']
halo = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1']
;,'1','1','1','1','1','1','1']
closerange = 1.0
farrange = 10.0
molecularH = 1
colorarr_ext = 0
mintform = 0.16503027; 0.17276607, 0.17792325
max_tform = 0.229737
compext = 'sb'

mfext = '.starball.mrk'

;-------------------------------------------------------------------

sfile = dir + file + '/' + file + '.starlog'
mfile = dir + file + '/steps/'+ file + '.' +steps[N_ELEMENTS(steps) - 1] + '.dir/' + file + '.' +steps[N_ELEMENTS(steps) - 1] ;+ '.halo.1' .starball.mrk'
tfiles = dir + file + '/steps/'+ file + '.' +steps    + '.dir/' + file + '.' +steps

;dir = '/astro/net/astro-quinn/fabio/REPOSITORY/e12Gals/h258.cosmo50cmb.1536gst1MbwKBH/'
;sfile = dir + 'h258.cosmo50cmb.1536gst1MbwKBH.starlog'
;mfile = dir + 'h258.cosmo50cmb.1536gst1MbwKBH.00456/h258.cosmo50cmb.1536gst1MbwKBH.00456'
;tfiles = dir + ['h258.cosmo50cmb.1536gst1MbwKBH.00455/h258.cosmo50cmb.1536gst1MbwKBH.00455','h258.cosmo50cmb.1536gst1MbwKBH.00456/h258.cosmo50cmb.1536gst1MbwKBH.00456']

rtipsy,mfile,header_mrk,g,d,s,/JUSTHEAD ;The header of the file the marked particles were taken from (probably the last output)
readarr,mfile+'.iord',header_mrk,iord_all_all,type = 'long',/ascii ;The marked particles
openr,1,mfile + mfext
readf,1,ntot,ngas,nstar
close,1
readcol,mfile+ mfext, mark,format = '(A)'
mrkind = LONG(mark[1:N_ELEMENTS(mark)-1]) - 1
read_tipsy_arr,mfile+'.iord',header_mrk,iord_all_all,type = 'long'
read_tipsy_arr,mfile+'.decomp',header_mrk,type = 'long'
;mass = track_mass(dir,file,steps,halo = halo,outplot = outplot)

FOR i = 0,N_ELEMENTS(steps) - 1 DO BEGIN
    tfile = tfiles[i]
    print,tfile
    IF KEYWORD_SET(outplot) THEN $
      track_starball, tfile, sfile, msol_per_sysmass, kpc_per_syslength, mrkind = mrkind, header_mrk = header_mrk, iord_all_all = iord_all_all, name = file + '.' +steps[i], halo = halo[i],farrange = farrange,closerange = closerange,colorarr_ext = colorarr_ext,mfext = mfext,molecularH = molecularH,outplot = compext, max_tform = max_tform, mintform = mintform ELSE $
      track_starball, tfile, sfile, msol_per_sysmass, kpc_per_syslength, mrkind = mrkind, header_mrk = header_mrk, iord_all_all = iord_all_all, name = file + '.' +steps[i], halo = halo[i],farrange = farrange,closerange = closerange,colorarr_ext = colorarr_ext,mfext = mfext,molecularH = molecularH, max_tform = max_tform, mintform = mintform 
ENDFOR
END

PRO track_starball, tfile, sfile, msol_per_sysmass, kpc_per_syslength_a, mrkind = mrkind, header_mrk = header_mrk, max_tform = max_tform, iord_all_all = iord_all_all, name = name, outplot = outplot, MOLECULARH = MOLECULARH,halo = halo,closerange = closerange,farrange = farrange,colorarr_ext = colorarr_ext,mfext = mfext,mintform = mintform;mfile_halo = mfile_halo,
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
IF NOT (keyword_set(mfext)) then mfext = '.halo.1.starball.mrk'

;--------------------------------------- Reading Data ---------------------------
rtipsy,tfile,header,g,d,s,/JUSTHEAD
;rtipsy,mfile,header_mrk,g,d,s,/JUSTHEAD
rtipsy,tfile+'.halo.'+halo,h,g,d,s 
print,'ntot: ',strtrim(h.n,2),', ndark ',strtrim(h.ndark,2),', ngas: ',strtrim(h.ngas,2),', nstar: ',strtrim(h.nstar,2)
IF KEYWORD_SET(mrkind) THEN BEGIN
;    readarr,tfile+'.halo.'+halo+'.iord',h,iord_gas,/ascii,part = 'gas',type = 'long'
;    readarr,tfile+'.halo.'+halo+'.iord',h,iord_star,/ascii,part = 'star',type = 'long'
;    IF KEYWORD_SET(colorarr_ext) THEN readarr,tfile+'.halo.'+halo+'.' + colorarr_ext,h,colorarr_all,/ascii,part = 'star',type = 'int'
;    readarr,mfile+'.iord',header_mrk,iord_all_all,type =
;    'long',/ascii
    read_tipsy_arr,tfile+'.halo.'+halo+'.iord',h,iord,type = 'long'
    iord_gas = iord[0:h.ngas - 1]
    iord_star = iord[h.ngas + h.ndark:h.n - 1]
    IF KEYWORD_SET(colorarr_ext) THEN BEGIN
        read_tipsy_arr,tfile+'.halo.'+halo+'.' + colorarr_ext,h,colorarr_all,type = 'int'
        colorarr_star = colorarr_all[h.ngas + h.ndark:h.n - 1]
    ENDIF
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

star_info = rstarlog(sfile,MOLECULARH = MOLECULARH)
IF NOT KEYWORD_SET(max_tform) THEN max_tform = MAX(star_info.timeform)
IF KEYWORD_SET(mintform) THEN BEGIN
    temp = MIN(ABS(star_info.timeform - mintform),minind)
    miniord = star_info[minind].IORDERSTAR
ENDIF
star_info = star_info[where(star_info.iorderstar ne 0)]
print,'Data Read'

;--------------------------------------- Selecting Stars and Gas in
;                                        Marked File --------------------------
IF KEYWORD_SET(mrkind) THEN BEGIN
    siord_all = iord_all_all[mrkind[where(mrkind gt header_mrk.ndark + header_mrk.ngas AND mrkind le header_mrk.ndark + header_mrk.ngas + header_mrk.nstar)]]

;Find which elements in the starlog correspond to the stars in the starball
;Select and array of starlog entries for the starball
;    sl_starball_ind = [0]
;    FOR is = 0L, N_ELEMENTS(siord_all) - 1 DO BEGIN
;        ind = (where(star_info.iorderstar eq siord_all[is]))[0]
;        sl_starball_ind = [sl_starball_ind,ind]
;    ENDFOR
;    sl_starball = star_info[sl_starball_ind[1:N_ELEMENTS(sl_starball_ind) - 1]]
;    match,siord_all,star_info.iorderstar,temp,sl_starball_ind
;    sl_starball = star_info[sl_starball_ind]

IF KEYWORD_SET(mintform) THEN siord_all = siord_all[where(siord_all ge miniord)]

;select the stars in the marked file that have already formed
    IF (MIN(siord_all) LT MAX(iord_star)) THEN BEGIN
         siord = siord_all[where(siord_all le max(iord_star))]
;        sind = fltarr(N_ELEMENTS(siord))
;        FOR is = 0, N_ELEMENTS(sind) - 1 DO BEGIN
;            sind[is] = where(siord_all[is] eq iord_star)
;        ENDFOR
;        s_starball = s[sind]
        match,siord,iord_star,temp,sind
;        match,siord,iord,temp,sind
        s_starball = s[sind]
        IF KEYWORD_SET(colorarr_ext) THEN colorarr = colorarr_star[sind]
        nostar = 0
    ENDIF ELSE BEGIN
        nostar = 1
    ENDELSE
    print,'# of star marked: ',STRTRIM(N_ELEMENTS(siord_all),2),', # of stars formed so far: ',STRTRIM(N_ELEMENTS(siord),2)

;select the gas in the marked file that stays gas
;    giord = iord_all_all[mrkind[where(mrkind le ngas)]]
    IF ((where(mrkind le header_mrk.ngas))[0] ne -1) THEN BEGIN
;        gind = fltarr(N_ELEMENTS(giord))
;        FOR ig = 0L, N_ELEMENTS(giord) - 1 DO BEGIN
;            IF (where(giord[ig] eq iord_gas))[0] NE -1 THEN gind[ig] = (where(giord[ig] eq iord_gas))[0]
;        ENDFOR
        match,giord,iord_gas,temp,gind
        g_starball= g[gind]        
        g_timeform = fltarr(N_ELEMENTS[g_starball]) 
        nogas = 0
    ENDIF ELSE BEGIN
        g_starball = -1
        g_timeform = -1
        nogas = 1
    ENDELSE

;select the gas in the starball that will become stars
    IF ((where(siord_all gt max(iord_star)))[0] ne -1) THEN BEGIN
        siord_future = siord_all[where(siord_all gt max(iord_star))]
        match,siord_future,star_info.iorderstar,temp,sl_future_ind ;match the stars in the starball that will form to the starlog
        sl_future = star_info[sl_future_ind]
        sfiord = sl_future.iordergas
        sftform = sl_future.timeform
        match,sfiord,iord_gas,sfind,gind   ;match the starlog starball information to the gas they formed from 
        gprior_starball = g[gind]
        gprior_timeform = sftform[sfind]
        IF g_starball[0] ne -1 THEN BEGIN ;Add gas that will form stars onto gas array
            g_starball = [g_starball,gprior_starball]
            g_timeform = [g_timeform,gprior_timeform]
            nogas = 0
        ENDIF ELSE BEGIN
            g_starball = gprior_starball
            g_timeform = gprior_timeform
        ENDELSE
    ENDIF

    IF 0 THEN BEGIN
        gprior_starball = [0]
        g_timeform = [0]
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
ENDIF

;----------------------------------------- Selecting Halo and Centering -------------
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

;---------------------------------------- Setting colors -------------------------
IF KEYWORD_SET(mrkind) THEN BEGIN
    IF KEYWORD_SET(colorarr_ext) THEN BEGIN
        mform = MAX(s.mass)
        ct = 39
        gas = 244
        disk = 50 ;Blue
        thdisk = 90 ; Light Blue
        pbbulge = 160 ; Green
        bulge = 180 ; Yellow
        halo_s = 215 ; Orange
        extra = 20 ; Purple
        IF NOT nostar THEN BEGIN
            s_color = colorarr
            IF ((where(colorarr eq 1))[0] ne -1 ) THEN s_color[where(colorarr eq 1)] = disk
            print,'Disk:         ',N_ELEMENTS([where(colorarr eq 1)]),N_ELEMENTS([where(colorarr eq 1)])/double(N_ELEMENTS(siord_all))
            IF ((where(colorarr eq 2))[0] ne -1 ) THEN s_color[where(colorarr eq 2)] = halo_s
            print,'Halo:         ',N_ELEMENTS([where(colorarr eq 2)]),N_ELEMENTS([where(colorarr eq 2)])/double(N_ELEMENTS(siord_all))
            IF ((where(colorarr eq 3))[0] ne -1 ) THEN s_color[where(colorarr eq 3)] = bulge
            print,'Bulge:        ',N_ELEMENTS([where(colorarr eq 3)]),N_ELEMENTS([where(colorarr eq 3)])/double(N_ELEMENTS(siord_all))
            IF ((where(colorarr eq 4))[0] ne -1 ) THEN s_color[where(colorarr eq 4)] = thdisk
            print,'Thick Disk:   ',N_ELEMENTS([where(colorarr eq 4)]),N_ELEMENTS([where(colorarr eq 4)])/double(N_ELEMENTS(siord_all))
            IF ((where(colorarr eq 5))[0] ne -1 ) THEN s_color[where(colorarr eq 5)] = pbbulge
            print,'Pseudob-bulge:',N_ELEMENTS([where(colorarr eq 5)]),N_ELEMENTS([where(colorarr eq 5)])/double(N_ELEMENTS(siord_all))
            IF ((where(colorarr eq 0))[0] ne -1 ) THEN s_color[where(colorarr eq 0)] = extra 
            print,'Extra:        ',N_ELEMENTS([where(colorarr eq 0)]),N_ELEMENTS([where(colorarr eq 0)])/double(N_ELEMENTS(siord_all))
        ENDIF
        g_color = intarr(N_ELEMENTS(g_starball)) + gas
        print,'Gas:          ',N_ELEMENTS(g_color),1.0 - N_ELEMENTS(colorarr)/double(N_ELEMENTS(siord_all))
        if N_ELEMENTS(g_starball) gt 1 THEN nogas = 0
    ENDIF ELSE BEGIN
        ct = 18
;g_color = FIX(g_timeform/max_tform*220.0) + 24
;g_color =  FIX((g_timeform - h.time)*170.0)
        IF KEYWORD_SET(mintform) THEN mintcolor = mintform ELSE mintcolor = MIN(g_timeform)
;        IF KEYWORD_SET(maxtform) THEN maxtcolor = maxtform ELSE maxtcolor = max_form
        g_color =  FIX((g_timeform - MIN(g_timeform))/(max_tform - MIN(mintcolor))*129.0) ;+ 124
        IF ((where(g_timeform eq 0))[0] ne -1 ) THEN BEGIN
            ind = where(g_timeform eq 0)
            g_color[ind] = 129
        ENDIF
        IF ((where(g_color gt 148.0))[0] ne -1 ) THEN BEGIN
            ind = where(g_color gt 129)
            g_color[ind] = 129
        ENDIF
;s_color = FIX(s_starball.tform/max_tform*220.0) + 24
;s_color = 299 - FIX(s_starball.tform/max_tform*170.0)
        IF NOT nostar THEN s_color = 130 + FIX((MAX(s_starball.tform) - s_starball.tform)/max_tform*130.0)
        IF N_ELEMENTS(g_starball) GT 1 THEN nogas = 0
    ENDELSE
ENDIF ELSE BEGIN
    s_color = 130 + FIX((MAX(s_halo.tform) - s_halo.tform)/max_tform*299.0)
ENDELSE

;---------------------------------------- Frame -----------------------------------
;!p.multi = [0,2,2]
!X.MARGIN = [6,1]
!Y.MARGIN = [3,2]
loadct,0
IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename = tfile+'.halo.1_image'+ outplot + '.eps',/times,ysize=18,xsize=18 ELSE window,0,xsize = 800,ysize = 800 ;392
multiplot,[2,2],gap=0.04,mtitle = name,/square,/doxaxis,/doyaxis
dxy = closerange           ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
plot,s_halo.x,s_halo.y,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],ytitle = 'Y [kpc]',/nodata;,title = name
;oplot,g_halo.x,g_halo.y,psym = 3,color = 150
oplot,s_halo.x,s_halo.y,psym = 3
loadct,ct
IF KEYWORD_SET(mrkind) THEN BEGIN
;oplot,g_starball.x,g_starball.y,psym = 3,color = 244
    IF NOT nogas THEN FOR ig = 0L,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].y,g_starball[ig].y],color = g_color[ig],psym = symcat(16),symsize = 0.2 ;psym = 3
    IF NOT nostar THEN FOR is = 0L,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].y,s_starball[is].y],color = s_color[is],psym = symcat(16),symsize = 0.2 ;psym = 3
ENDIF ELSE FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].y,s_halo[is].y],psym = 3,color = s_color[is]    
multiplot,/doxaxis,/doyaxis

dxy = farrange          ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
loadct,0
plot,s_halo.x,s_halo.y,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],/nodata;,title = name
oplot,g_halo.x,g_halo.y,psym = 3,color = 150
oplot,s_halo.x,s_halo.y,psym = 3
loadct,ct
IF KEYWORD_SET(mrkind) THEN BEGIN
    IF NOT nogas THEN FOR ig = 0L,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].y,g_starball[ig].y],psym = 3,color = g_color[ig]
    IF NOT nostar THEN FOR is = 0L,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].y,s_starball[is].y],psym = 3,color = s_color[is]
ENDIF ELSE  FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].y,s_halo[is].y],psym = 3,color = s_color[is]
multiplot,/doxaxis,/doyaxis    


;dxy = closerange           ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
;loadct,0
;plot,s_halo.x,s_halo.z,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],xtitle = 'X [kpc]',ytitle = 'Z [kpc]',/nodata;,title = name
;oplot,g_halo.x,g_halo.z,psym = 3,color = 150
;oplot,s_halo.x,s_halo.z,psym = 3
;loadct,ct
;IF KEYWORD_SET(mrkind) THEN BEGIN
;    IF NOT nogas THEN FOR ig = 0L,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].z,g_starball[ig].z],psym = 3,color = g_color[ig]
;    IF NOT nostar THEN FOR is = 0L,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].z,s_starball[is].z],psym = 3,color = s_color[is]
;ENDIF ELSE FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].z,s_halo[is].z],psym = 3,color = s_color[is]
;multiplot,/doxaxis,/doyaxis


dxy = closerange           ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
loadct,0
plot,s_halo.x,s_halo.z,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],xtitle = 'X [kpc]',ytitle = 'Y [kpc]',/nodata;,title = name
oplot,g_halo.x,g_halo.y,psym = 3;,color = 150
;oplot,s_halo.x,s_halo.z,psym = 3
loadct,ct
IF KEYWORD_SET(mrkind) THEN BEGIN
    IF NOT nogas THEN FOR ig = 0L,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].y,g_starball[ig].y],color = g_color[ig],psym = symcat(16),symsize = 0.2 ;psym = 3
    IF NOT nostar THEN FOR is = 0L,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].y,s_starball[is].y],color = s_color[is],psym = symcat(16),symsize = 0.2 ;psym = 3
ENDIF ELSE FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].z,s_halo[is].z],psym = 3,color = s_color[is]
multiplot,/doxaxis,/doyaxis
    
dxy = farrange          ;100; 10.0*MEAN([STDDEV(s_halo.x),STDDEV(s_halo.y)])
loadct,0
plot,s_halo.x,s_halo.z,psym = 3,xrange = [-1.0*dxy, dxy],yrange = [-1.0*dxy, dxy],xtitle = 'X [kpc]',/nodata;,title = name
oplot,g_halo.x,g_halo.z,psym = 3,color = 150
oplot,s_halo.x,s_halo.z,psym = 3
loadct,ct
IF KEYWORD_SET(mrkind) THEN BEGIN
    IF NOT nogas THEN FOR ig = 0L,N_ELEMENTS(g_starball)-1 DO oplot,[g_starball[ig].x,g_starball[ig].x],[g_starball[ig].z,g_starball[ig].z],psym = 3,color = g_color[ig]
    IF NOT nostar THEN FOR is = 0L,N_ELEMENTS(s_starball)-1 DO oplot,[s_starball[is].x,s_starball[is].x],[s_starball[is].z,s_starball[is].z],psym = 3,color = s_color[is]
ENDIF ELSE FOR is = 0L,N_ELEMENTS(s_halo)-1 DO oplot,[s_halo[is].x,s_halo[is].x],[s_halo[is].z,s_halo[is].z],psym = 3,color = s_color[is]

IF KEYWORD_SET(outplot) THEN  device,/close ;else stop
loadct,0
multiplot,/reset

;---------------------- SFH ------------------------------------------------------------------------
loadct,39
!p.multi = 0
!X.MARGIN = [12,3]
!Y.MARGIN = [6,2]
IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename=tfile+'.halo.1_sfh.eps',/times,ysize=18,xsize=18 ELSE window,2,xsize = 630,ysize = 392
sfr,s_halo,massunit = msol_per_sysmass,timeunit = timeunit,binsize = 1e7,xrange = [6,9],yrange = [0,0.25],xrnage = [0,timeunit*max_tform/1e9],massform = MAX(s.mass),thick = 3
;sfr,s_halo,massunit = msol_per_sysmass,timeunit = timeunit,binsize = 1e7,xrange = [6,8],yrange = [0,0.05],xrnage = [0,timeunit*max_tform/1e9],massform = MAX(s.mass)
IF (KEYWORD_SET(mrkind)) THEN $
  IF NOT nostar THEN sfr,s_starball,massunit = msol_per_sysmass,timeunit = timeunit,/overplot,color = 240,binsize = 1e7,thick = 3
IF KEYWORD_SET(outplot) THEN device,/close else stop

;---------------------- Energy  ------------------------------------------
IF 0 then begin
d_halo.x = d_halo.x/kpc_per_syslength
d_halo.y = d_halo.y/kpc_per_syslength
d_halo.z = d_halo.z/kpc_per_syslength
g_halo.x = g_halo.x/kpc_per_syslength
g_halo.y = g_halo.y/kpc_per_syslength
g_halo.z = g_halo.z/kpc_per_syslength
s_halo.x = s_halo.x/kpc_per_syslength
s_halo.y = s_halo.y/kpc_per_syslength
s_halo.z = s_halo.z/kpc_per_syslength
s_starball.x = s_starball.x/kpc_per_syslength
s_starball.y = s_starball.y/kpc_per_syslength
s_starball.z = s_starball.z/kpc_per_syslength

Rxyz = sqrt(s_halo.x*s_halo.x+s_halo.y*s_halo.y+s_halo.z*s_halo.z) 
v2=(s_halo.vx*s_halo.vx+s_halo.vy*s_halo.vy+s_halo.vz*s_halo.vz)
ke=0.5*v2
Etot=(ke+s_halo.phi)
max_etot = MAX(Etot)
s_halo.phi=s_halo.phi-max_etot
Etot=(ke+s_halo.phi)

s_starball_save = s_starball
v2_starball=(s_starball.vx*s_starball.vx+s_starball.vy*s_starball.vy+s_starball.vz*s_starball.vz)
ke_starball = 0.5*v2_starball
s_starball.phi = s_starball.phi - max_etot
Etot_starball = (ke_starball + s_starball.phi)
    
Jz = (s_halo.x*s_halo.vy-s_halo.y*s_halo.vx)
Jz_starball = (s_starball.x*s_starball.vy-s_starball.y*s_starball.vx)
Jx = (s_halo.y*s_halo.vz-s_halo.z*s_halo.vy)
Jy = (s_halo.z*s_halo.vx-s_halo.x*s_halo.vz)
J=sqrt(Jx*Jx+Jy*Jy+Jz*Jz)

;rotfile = 'h258.cosmo50cmb.1536gst1MbwKBH.00455.rot.1'
;lines=file_lines(rotfile)
;readcol,rotfile,radius,skipline=lines-2,/silent
;radius=radius[0]
;nstars = n_elements(s_halo.x)
;JzJzmax=    replicate(-2.,nstars)

;getlcirc,'dummy',h,g_halo,d_halo,s_halo,E,Jc,rmax=radius,zmax=0.2/kpc_per_syslength,rotfile=rotfile
;Energy = Etot
;sort=sort(Energy)
;Energy = Energy[sort]
;Jzmax = spline(E,Jc,Energy)
;sort=sort(sort)
;Jzmax = Jzmax[sort]
;Energy = Energy[sort]
;keep = where(Jzmax lt max(Jc) and Energy lt max(E)) ; get rid of weird outputs
;Jzmax=Jzmax[keep]
;Energy=Energy[keep]
;region=region[keep]
;JzJzmax=Jz/Jzmax
;;disk=region[where(JzJzmax[region] gt .8 and JzJzmax[region] lt 1.1)]
ENDIF

IF 0 THEN BEGIN
IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename = tfile+'.energy'+ outplot + '.eps',/times,ysize=18,xsize=18 ELSE window,1,xsize = 800,ysize = 800 ;392
multiplot,[1,2],gap=0.04,mtitle = name,/doxaxis,/doyaxis
plot,Etot,Jz,psym=3,title='Energy',xtitle='Etot',ytitle='J_z' ; check the range for Ecut
FOR is = 0L,N_ELEMENTS(s_starball)-1 DO  oplot,[Etot_starball[is],Etot_starball[is]], [Jz_starball[is],Jz_starball[is]],color = s_color[is],psym = 3

multiplot
plot,Etot,Jz,psym=3,title='Energy',xtitle='Etot',ytitle='J_z',yrange = [-0.00001,0.0001],xrange = [-0.22,-0.03] ; check the range for Ecut
FOR is = 0L,N_ELEMENTS(s_starball)-1 DO  oplot,[Etot_starball[is],Etot_starball[is]], [Jz_starball[is],Jz_starball[is]],color = s_color[is],psym = 3
ENDIF
multiplot,/reset
close,/all
END
