;This reads in a starlog and plots the star formation history  for
;different sections
PRO track_sfh_master
dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/'
file = 'h799.cosmo25cmb.3072g1MBWK'
steps =['00024','00036','00048','00060','00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00512']
halo = ['1',    '1',    '1',    '1',    '1',    '2',    '1',    '2',    '2',    '1',    '1',    '1',    '1',    '1',    '1',    '1',    '1']

mfile = dir + file + '/steps/'+ file + '.' +steps[N_ELEMENTS(steps) - 1] + '.dir/' + file + '.' +steps[N_ELEMENTS(steps) - 1] ;+ '.halo.' + halo[N_ELEMENTS(steps) - 1]
rtipsy,mfile,h,g,d,s,/JUSTHEAD  ;Read information about the halo
readarr,mfile+'.iord',h,iord,type = 'long',/ascii ;read the iords for the halo

;component 0 (selected from a grp file, i.e. a halo)
readarr,mfile+'.amiga.grp',h,ind_c0,/ascii,type = 'star'
siord_c0 = iord[where(ind_c0 eq FIX(halo[N_ELEMENTS(steps) - 1]))]

;component 1 (a marked file)
openr,1,mfile+'.halo.1.starball.mrk'
readf,1,ntot,ngas,nstar
close,1
readcol,mfile+'.halo.1.starball.mrk',mark ;read the marked particles
ind_c1 = mark[1:N_ELEMENTS(mark)-1] - 1 ;Adjust for +1 index
siord_c1 = iord[ind_c1[where(ind_c1 gt h.ndark + h.ngas AND ind_c1 le h.ndark + h.ngas + h.nstar)]] 
readcol,dir + file + '/' + file + '.log',dTime, z, E, T, U, Eth, Lx, Ly, Lz, WallTime, dWMax, dImax, dEMax, dMultiEff
units = tipsyunits(dir + file + '/' + file + '.param')
time = fltarr(N_ELEMENTS(steps))
FOR i = 0,N_ELEMENTS(steps) - 1 DO BEGIN
    tfile = dir + file + '/steps/'+ file + '.' +steps[i] + '.dir/' + file + '.' +steps[i]
    rtipsy,tfile,h,g,d,s,/justhead
    temp = min(abs(z - (1.0/h.time - 1.0)),ind)
    time[i] = dTime[ind]*units.timeunit
ENDFOR
outplot =  dir + file + '/steps/'+ file + '.' +steps + '.dir/' + file + '.' +steps + '.halo.1_sfh.eps'
track_sfh, dir, file, units, siord_c0, comp1 = siord_c1, max_tform = time, outplot = outplot

END

PRO track_sfh, dir, file, units, comp0, comp1 = comp1, comp2 = comp2, outplot = outplot, molecularH = molecularH, max_tform = max_tform
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
    thick = 2
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    thick = 1
ENDELSE

sfile = dir + file + '/' + file + '.starlog'
star_info = rstarlog(sfile, MOLECULARH = MOLECULARH)
IF NOT KEYWORD_SET(max_tform) THEN max_tform = MAX(star_info.timeform)
star_info = star_info[where(star_info.iorderstar ne 0)]

;sl_c0_ind = [0]
;FOR is = 0L, N_ELEMENTS(comp0) - 1 DO BEGIN
;    ind = (where(star_info.iorderstar eq sl_c0_ind[is]))[0]
;    sl_c0_ind = [sl_c0_ind,ind]
;ENDFOR
;sl_c0 = star_info[sl_c0_ind[1:N_ELEMENTS(sl_c0_ind) - 1]]
match,comp0,star_info.iorderstar,temp,sl_c0_ind
sl_c0 = star_info[sl_c0_ind]

;Component 1 ------------------------------------
IF KEYWORD_SET(comp1) THEN BEGIN
;    sl_c1_ind = [0]
;    FOR is = 0L, N_ELEMENTS(comp1) - 1 DO BEGIN
;        ind = (where(star_info.iorderstar eq sl_c1_ind[is]))[0]
;        sl_c1_ind = [sl_c1_ind,ind]
;    ENDFOR
;    sl_c1 = star_info[sl_c1_ind[1:N_ELEMENTS(sl_c1_ind) - 1]]
    match,comp1,star_info.iorderstar,temp,sl_c1_ind
    sl_c1 = star_info[sl_c1_ind]
    color1 = 240
ENDIF

;Component 2 ------------------------------------
IF KEYWORD_SET(comp2) THEN BEGIN
;    sl_c2_ind = [0]
;    FOR is = 0L, N_ELEMENTS(comp2) - 1 DO BEGIN
;        ind = (where(star_info.iorderstar eq sl_c2_ind[is]))[0]
;        sl_c2_ind = [sl_c2_ind,ind]
;    ENDFOR
;    sl_c2 = star_info[sl_c2_ind[1:N_ELEMENTS(sl_c2_ind) - 1]]
    match,comp2,star_info.iorderstar,temp,sl_c2_ind
    sl_c2 = star_info[sl_c2_ind]
    color2 = 50
ENDIF

FOR i = 0, N_ELEMENTS(max_tform) - 1 DO BEGIN
    IF KEYWORD_SET(outplot) THEN device,/color,bits_per_pixel=8,filename = outplot[i],/times,ysize=12,xsize=18 ELSE window,1
    sl_c0_cut = sl_c0[where(sl_c0.timeform*units.timeunit le max_tform[i])]
    IF KEYWORD_SET(comp1) THEN sl_c1_cut = sl_c1[where(sl_c1.timeform*units.timeunit le max_tform[i])]
    IF KEYWORD_SET(comp2) THEN sl_c1_cut = sl_c1[where(sl_c1.timeform*units.timeunit le max_tform[i])]
    sfr,sl_c0_cut,/starlog, massform= units.ISTARMASS,timeunit = units.timeunit,xrange = [0,MAX(star_info.timeform)*units.timeunit/1e9],massunit = 1,binsize = 1e8,yrange = [0,0.06],thick = thick
    IF KEYWORD_SET(comp1) THEN sfr,sl_c1_cut,/starlog,color = color1,massform = units.ISTARMASS,timeunit = units.timeunit,massunit = 1,binsize = 1e8,thick = thick,/overplot
    IF KEYWORD_SET(comp2) THEN sfr,sl_c2_cut,/starlog,color = color2,massform = units.ISTARMASS,timeunit = units.timeunit,massunit = 1,binsize = 1e8,thick = thick,/overplot
    IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
ENDFOR
END
