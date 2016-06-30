;dir = '/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_1e5/H2.repo.7.20.old'
;pfile = '/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_1e5/H2.repo.7.20.old/MW_disk.H2.param'
;step = '00010'
;filename = 'Disk_Collapse_1e5.00150_H2'

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/'
;filename = 'h603.cosmo50cmb.3072g14HBWK'
;step = '00408'
;pfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param'

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/'
;pfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param'
;step = '00408'
;step = '00504'
;step = '00276'
;step = '00288'
;filename = 'h986.cosmo50cmb.3072g14HBWK'

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/'
;pfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param'
;step = '00504'
;step = '00288'
;filename = 'h986.cosmo50cmb.3072gs1MbwK'
;title = 'Metals, no ' + textoidl('H_2')

;dr = 10.0
;maxr = 10.0

PRO redistribution_prep,dir,filename,pfile,step,dr = dr,maxr = maxr, maxt = maxt,title = title
IF NOT KEYWORD_SET(dr) THEN dr = 2.0
IF NOT KEYWORD_SET(maxr) THEN maxr = 8.0
IF NOT KEYWORD_SET(maxt) THEN maxt = 1.2e4
cd,dir
filebase = filename + '.' + step + '.dir/' + filename + '.' + step + '.halo.1'
;filebase = filename + '.' + step
rtipsy,filebase,h,g,d,s
readarr,filebase + '.iord',h,iord,/ascii,part = 'gas'
units = tipsyunits(pfile)
fixunits,h,g,d,s,units
grad = SQRT(g.x*g.x + g.y*g.y)

ind = where(g.tempg LE maxt AND ABS(g.z) LT 4) 
iord_set = iord[ind]
writecol,filename + '.' + step + '.cold.iord',LONG(iord_set)

FOR i = 0, maxr/dr DO BEGIN
    ind = where(grad GE i*dr AND grad LT (i + 1)*dr AND g.tempg LE maxt AND ABS(g.z) LT 4) 
    iord_set = iord[ind]
    writecol,filename + '.' + step + '.' + STRTRIM(FIX((i+1)*dr),2) + '.iord',LONG(iord_set)
ENDFOR



END

;sfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.starlog'
;steps = ['00408','00512']

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/'
;filename = 'h986.cosmo50cmb.3072g14HBWK'
;pfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.param'
;sfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/h986.cosmo50cmb.3072g14HBWK.starlog'
;steps = ['00408','00512']
;steps = ['00504','00512']
;steps = ['00276','00288']
;steps = ['00288','00300']
;center0 = [0,0,0] ;276
;center1 = [0,0,0] ;288
;center2 = [0,0,0] ;300
;title = textoidl('H_2')

;dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/steps/'
;filename = 'h986.cosmo50cmb.3072gs1MbwK'
;pfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.param'
;sfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.starlog.
;steps = ['00408','00512']
;steps = ['00504','00512']
;steps = ['00276','00288']
;steps = ['00288','00300']
;center0 = [-9.656405e-6,-5.117131e-6,7.207145e-7] ;276
;center1 = [0,0,0] ;288
;center2 = [0,0,0]
;title = 'Metals, no ' + textoidl('H_2')

;steps = ['00010','00020']
;sfile = '/astro/store/student-scratch1/christensen/MolecH/12M/Disk_Collapse_1e5/H2.repo.7.20.old/Disk_Collapse_1e5.00150_H2.starlog'

;steps = ['00504','00512']

;maxr = 30
;dr = 2
PRO redistribution,dir,filename,pfile,sfile,steps,center0 = center0, center1 = center1, dr = dr,maxr = maxr, maxt = maxt,MOLECULARH = MOLECULARH,nostarlog = nostarlog,title = title
set_plot,'x'
loadct,39
formatplot
IF NOT KEYWORD_SET(dr) THEN dr = 2.0
IF NOT KEYWORD_SET(maxr) THEN maxr = 8.0
IF NOT KEYWORD_SET(maxt) THEN maxt = 1.2e4
color = (findgen(maxr/dr + 1) + 1)*240/(maxr/dr + 1)

cd,dir
filebase0 = filename + '.' + steps[0] + '.dir/' + filename + '.' + steps[0] + '.halo.1'
filebase1 = filename + '.' + steps[1] + '.dir/' + filename + '.' + steps[1] + '.halo.1'
;filebase0 = filename + '.' + steps[0]
;filebase1 = filename + '.' + steps[1] 
rtipsy,filebase0,h0,g0,d0,s0
rtipsy,filebase1,h1,g1,d1,s1
g0.x = g0.x - center0[0]
g0.y = g0.y - center0[1]
g0.z = g0.z - center0[2]
s0.x = s0.x - center0[0]
s0.y = s0.y - center0[1]
s0.z = s0.z - center0[2]
d0.x = d0.x - center0[0]
d0.y = d0.y - center0[1]
d0.z = d0.z - center0[2]

g1.x = g1.x - center1[0]
g1.y = g1.y - center1[1]
g1.z = g1.z - center1[2]
s1.x = s1.x - center1[0]
s1.y = s1.y - center1[1]
s1.z = s1.z - center1[2]
d1.x = d1.x - center1[0]
d1.y = d1.y - center1[1]
d1.z = d1.z - center1[2]

units = tipsyunits(pfile)
fixunits,h0,g0,d0,s0,units
fixunits,h1,g1,d1,s1,units
readarr,filebase0 + '.iord',h0,iord0_gas,/ascii,part = 'gas'
readarr,filebase1 + '.iord',h1,iord1,/ascii
;eadarr,filebase1 + '.coolontime',h1,coolontime,/ascii,part = 'gas'
iord1_gas  = iord1[0:h1.ngas - 1]
iord1_star = iord1[h1.ngas + h1.ndark: h1.n - 1]
IF NOT KEYWORD_SET(nostarlog) THEN star_info_all = rstarlog(sfile,MOLECULARH = MOLECULARH)
maxtform0 = MAX(s0.tform)
maxtform1 = MAX(s1.tform)
IF NOT KEYWORD_SET(nostarlog) THEN  star_info = star_info_all[where(star_info_all.timeform gt maxtform0 AND star_info_all.timeform le maxtform1)] ;select for stars that had not formed by step 0

;FOR i = 0, maxr/dr DO BEGIN
readcol,filename + '.' + steps[0] + '.cold.iord',iord_gas

;Match Gas
;    match,iord_gas,iord0_gas,tempg0,ind0_gas ;Find the indicies of the gas particles from step0 that have the same iord values of those in the file
match,iord_gas,iord1_gas,tempg1,ind1_gas ;Find the indicies of the gas particles from     step1 that have the same iord values of those in the file
gmatch1 = g1[ind1_gas] ;sort gas particles so that it is in the same order as iord_gas[tempg1]
;coolontime1 = coolontime[ind1_gas]
iord_gas_gas = iord_gas[tempg1] ;Iord values of gas particles that last the two steps

match,iord_gas_gas,iord0_gas,tempg0,ind0_gas
ggmatch0 = g0[ind0_gas]
gmatch1 = gmatch1[tempg0]
;coolontime1 = coolontime1[tempg0]

;Match Stars
IF NOT KEYWORD_SET(nostarlog) THEN BEGIN
    match,iord_gas,star_info.iordergas,temps1,ind1_sl ;Find the indicies of the starlog file data that represent stars from from gas with the same iord values of those in the file
    star_info_gas = star_info[ind1_sl] ;star log data for the gas that becomes stars
    iord_gas_star = iord_gas[temps1] ;iord values of the gas that turns into stars

    match,iord_gas_star,iord0_gas,tempg0,ind0_gas ; match ga particles to their iords
    iord_gas_star = iord_gas_star[tempg0] ;fix order to iords
    star_info_gas = star_info_gas[tempg0] ;fix order of star log info
    iord0_gas_star = iord0_gas[ind0_gas] ;fix order of selected gas iords
    gsmatch0 = g0[ind0_gas]     ;select gas particle

    match,star_info_gas.IORDERSTAR,iord1_star,temp,ind1_star ;match iords of formed stars to stellar iord
    smatch1 = s1[ind1_star]
    gsmatch0 = gsmatch0[temp]

    gmatch0 = [ggmatch0,gsmatch0]
ENDIF ELSE gmatch0 = ggmatch0

IF 1 THEN BEGIN
    angmom, ggmatch0, jvecgg0, lvec0, jmaggg0, ltot0
    angmom, gmatch1, jvecg1, lvec1, jmagg1, ltot1
    angmom, gmatch0, jvecg0, lvec0, jmagg0, ltot0
    IF NOT KEYWORD_SET(nostarlog) THEN angmom, gsmatch0, jvecgs0, lvec0, jmaggs0, ltot0
    IF NOT KEYWORD_SET(nostarlog) THEN angmom, smatch1, jvecs1, lvec1, jmags1, ltot1
    gg0rad = [[ggmatch0.x],[ggmatch0.y],[ggmatch0.z]]
    IF NOT KEYWORD_SET(nostarlog) THEN gs0rad = [gsmatch0.x,gsmatch0.y,gsmatch0.z]
    g0rad = [[gmatch0.x],[gmatch0.y],[gmatch0.z]]
    g1rad = [[gmatch1.x],[gmatch1.y],[gmatch1.z]]
    IF NOT KEYWORD_SET(nostarlog) THEN s1rad = [[smatch1.x],[smatch1.y],[smatch1.z]]
    gg0vel = [[ggmatch0.vx],[ggmatch0.vy],[ggmatch0.vz]]
    IF NOT KEYWORD_SET(nostarlog) THEN gs0vel = [[gsmatch0.vx],[gsmatch0.vy],[gsmatch0.vz]]
    g0vel = [[gmatch0.vx],[gmatch0.vy],[gmatch0.vz]]
    g1vel = [[gmatch1.vx],[gmatch1.vy],[gmatch1.vz]]
    IF NOT KEYWORD_SET(nostarlog) THEN s1vel = [[smatch1.vx],[smatch1.vy],[smatch1.vz]]

;    ggr0 = SQRT(ggmatch0.x*ggmatch0.x + ggmatch0.y*ggmatch0.y)
;    IF NOT KEYWORD_SET(nostarlog)THEN gsr0 = SQRT(gsmatch0.x*gsmatch0.x + gsmatch0.y*gsmatch0.y)
;    gr0 = SQRT(gmatch0.x*gmatch0.x + gmatch0.y*gmatch0.y)
;    gr1 = SQRT(gmatch1.x*gmatch1.x + gmatch1.y*gmatch1.y)
;    IF NOT KEYWORD_SET(nostarlog) THEN gs1 = SQRT(smatch1.x*smatch1.x + smatch1.y*smatch1.y)
    ggr0 = SQRT(ggmatch0.x*ggmatch0.x + ggmatch0.y*ggmatch0.y + ggmatch0.z*ggmatch0.z)
    IF NOT KEYWORD_SET(nostarlog)THEN gsr0 = SQRT(gsmatch0.x*gsmatch0.x + gsmatch0.y*gsmatch0.y + gsmatch0.z*gsmatch0.z)
    gr0 = SQRT(gmatch0.x*gmatch0.x + gmatch0.y*gmatch0.y + gmatch0.z*gmatch0.z)
    gr1 = SQRT(gmatch1.x*gmatch1.x + gmatch1.y*gmatch1.y + gmatch1.z*gmatch1.z)
    IF NOT KEYWORD_SET(nostarlog) THEN gsr1 = SQRT(smatch1.x*smatch1.x + smatch1.y*smatch1.y + smatch1.z*smatch1.z)
   
    ggvr0 = REFORM(TRANSPOSE(dot(TRANSPOSE(gg0rad),TRANSPOSE(gg0vel))/mag(gg0rad)))
    IF NOT KEYWORD_SET(nostarlog)THEN gsvr0 =  REFORM(TRANSPOSE(dot(TRANSPOSE(gs0rad),TRANSPOSE(gs0vel))/mag(gs0rad)))
    gvr0 =  REFORM(TRANSPOSE(dot(TRANSPOSE(g0rad),TRANSPOSE(g0vel))/mag(g0rad)))
    gvr1 =  REFORM(TRANSPOSE(dot(TRANSPOSE(g1rad),TRANSPOSE(g1vel))/mag(g1rad)))
    IF NOT KEYWORD_SET(nostarlog) THEN gsr1 =  REFORM(TRANSPOSE(dot(TRANSPOSE(s1rad),TRANSPOSE(s1vel))/mag(s1rad)))

    hdr = 0.5
    hbins = findgen(maxr/hdr)*hdr + hdr/2.0
    hbinst = fltarr(maxr/hdr)

    hbinsl = fltarr(maxr/hdr)
    hbinslh = fltarr(maxr/hdr)
    hbinslc = fltarr(maxr/hdr)
    hbinsls = fltarr(maxr/hdr)
    hbinslsdv = fltarr(maxr/hdr)
    hbinslhsdv = fltarr(maxr/hdr)
    hbinslcsdv = fltarr(maxr/hdr)
    hbinslssdv = fltarr(maxr/hdr)

    hbinsr = fltarr(maxr/hdr)
    hbinsrh = fltarr(maxr/hdr)
    hbinsrc = fltarr(maxr/hdr)
    hbinsrs = fltarr(maxr/hdr)
    hbinsrsdv = fltarr(maxr/hdr)
    hbinsrhsdv = fltarr(maxr/hdr)
    hbinsrcsdv = fltarr(maxr/hdr)
    hbinsrssdv = fltarr(maxr/hdr)

    hbinsv = fltarr(maxr/hdr)
    hbinsvh = fltarr(maxr/hdr)
    hbinsvc = fltarr(maxr/hdr)
    hbinsvs = fltarr(maxr/hdr)
    hbinsvsdv = fltarr(maxr/hdr)
    hbinsvhsdv = fltarr(maxr/hdr)
    hbinsvcsdv = fltarr(maxr/hdr)
    hbinsvssdv = fltarr(maxr/hdr)

    colorsr = (((findgen(maxr/hdr) + 1)/(maxr/hdr)))*240
    window,0,xsize = 600,ysize = 600
;   IF N_ELEMENTS(indr) gt 1 THEN 
    plot,gmatch1.x,gmatch1.y,psym = 3,xrange = [-30,30],yrange = [-30,30] ;ELSE $
;     plot,[gmatch1.x,gmatch1.x],[gmatch1.y,gmatch1.y],psym = 3,xrange = [-30,30],yrange = [-30,30]   

;,xrange = [-20,20],yrange = [-20,20] 
;    plot,ggmatch0.x,ggmatch0.y,psym = 3
;    stop
    FOR i = 0, maxr/hdr - 1 DO BEGIN
        indr = where(ggr0 ge i*hdr AND $
                     ggr0 lt (i+1)*hdr); AND $
;                     jvecgg0[*,2] ge 0 and jvecg1[*,2] ge 0)
        IF NOT KEYWORD_SET(nostarlog) THEN indr_s = where(gsr0 ge i*hdr AND $
                     gsr0 lt (i+1)*hdr); AND $
;                     jvecgg0[*,2] ge 0 and jvecg1[*,2] ge 0)
        indr_hot = where(ggr0 ge i*hdr AND $
                         ggr0 lt (i+1)*hdr AND $
                         (gmatch1.tempg ge maxt)); OR $
;                         coolontime gt 0)); AND $
;                         jvecgg0[*,2] ge 0 and jvecg1[*,2] ge 0)
        indr_cold = where(ggr0 ge i*hdr AND $
                          ggr0 lt (i+1)*hdr AND $
                          (gmatch1.tempg lt maxt)); AND $
;                          coolontime le 0)) ;AND $
;                          jvecgg0[*,2] ge 0 and jvecg1[*,2] ge 0)
;        stop
        IF N_ELEMENTS(indr) gt 1 THEN oplot,gmatch1[indr].x,gmatch1[indr].y,psym = 3,color = colorsr[i] ;LSE oplot,[gmatch1[indr].x,gmatch1[indr].x],[gmatch1[indr].y,gmatch1[indr].y],psym = 3,color = colorsr[i]
        IF (indr_hot)[0] ne -1 THEN BEGIN
            hbinst[i] = TOTAL(gmatch1[indr_hot].mass)/TOTAL(gmatch1[indr].mass) 
            hbinslh[i] = MEDIAN([(jvecg1[indr_hot,2] - jvecgg0[indr_hot,2])/jvecgg0[indr_hot,2]])
            IF N_ELEMENTS(indr_hot) gt 1 THEN hbinslhsdv[i] = STDEV([(jvecg1[indr_hot,2] - jvecgg0[indr_hot,2])/jvecgg0[indr_hot,2]]) ELSE hbinsrhsdv[i] = 0
            hbinsrh[i] = MEDIAN(gr1[indr_hot])
            IF N_ELEMENTS(indr_hot) gt 1 THEN hbinsrhsdv[i] = STDEV([gr1[indr_hot]]) ELSE hbinsrhsdv[i] = 0
            hbinsvh[i] = MEDIAN(gvr1[indr_hot])
            IF N_ELEMENTS(indr_hot) gt 1 THEN hbinsvhsdv[i] = STDEV([gvr1[indr_hot]]) ELSE hbinsvhsdv[i] = 0

            IF N_ELEMENTS(indr_hot) gt 1 THEN oplot,gmatch1[indr_hot].x,gmatch1[indr_hot].y,psym = 3,color = colorsr[i] ELSE oplot,[gmatch1[indr_hot].x,gmatch1[indr_hot].x],[gmatch1[indr_hot].y,gmatch1[indr_hot].y],psym = 3,color = colorsr[i]
        ENDIF ELSE BEGIN
            hbinst[i] = 0
            hbinslh[i] = 0
            hbinsrh[i] = 0
            hbinsrhsdv[i] = 0
            hbinsvh[i] = 0
            hbinsvhsdv[i] = 0
        ENDELSE
        IF (indr_cold)[0] ne -1 THEN BEGIN
            hbinslc[i] = MEDIAN((jvecg1[indr_cold,2] - jvecgg0[indr_cold,2])/jvecgg0[indr_cold,2])
            IF N_ELEMENTS(indr_cold) gt 1 THEN hbinslcsdv[i] = STDEV(jvecg1[indr_cold,2] - jvecgg0[indr_cold,2]) ELSE hbinslcsdv[i] = 0
            hbinsrc[i] = MEDIAN(gr1[indr_cold])
            IF N_ELEMENTS(indr_cold) gt 1 THEN hbinsrcsdv[i] = STDEV(gr1[indr_cold])ELSE hbinsrcsdv[i] = 0
            hbinsvc[i] = MEDIAN(gvr1[indr_cold])
            IF N_ELEMENTS(indr_cold) gt 1 THEN hbinsvcsdv[i] = STDEV(gvr1[indr_cold])ELSE hbinsvcsdv[i] = 0
        ENDIF ELSE BEGIN
            hbinsl[i] = 0
            hbinslsdv[i] = 0
            hbinsr[i] = 0
            hbinsrsdv[i] = 0
            hbinsv[i] = 0
            hbinsvsdv[i] = 0
        ENDELSE
        IF (indr)[0] ne -1 THEN BEGIN
            hbinsl[i] = MEDIAN((jvecg1[indr,2] - jvecgg0[indr,2])/jvecgg0[indr,2])
            IF N_ELEMENTS(indr) gt 1 THEN hbinslsdv[i] = STDEV([jvecg1[indr,2] - jvecgg0[indr,2]])ELSE hbinslsdv[i] = 0
            hbinsr[i] = MEDIAN(gr1[indr])
            IF N_ELEMENTS(indr) gt 1 THEN hbinsrsdv[i] = STDEV([gr1[indr]])ELSE  hbinsrsdv[i] = 0
            hbinsv[i] = MEDIAN(gvr1[indr])
            IF N_ELEMENTS(indr) gt 1 THEN hbinsvsdv[i] = STDEV([gvr1[indr]])ELSE  hbinsvsdv[i] = 0
        ENDIF ELSE BEGIN
            hbinsl[i] = 0
            hbinslsdv[i] = 0
            hbinsr[i] = 0
            hbinsrsdv[i] = 0
            hbinsv[i] = 0
            hbinsvsdv[i] = 0
        ENDELSE
        IF NOT KEYWORD_SET(nostarlog) THEN BEGIN
        IF (indr_s)[0] ne -1 THEN BEGIN
            hbinsls[i] = MEDIAN((jvecs1[indr_s,2] - jvecgs0[indr_s,2])/jvecgs0[indr_s,2])
            IF N_ELEMENTS(indr_s) gt 1 THEN hbinslssdv[i] = STDEV([jvecs1[indr_s,2] - jvecgs0[indr_s,2]]) ELSE hbinslssdv[i] = 0
            hbinsrs[i] = MEDIAN(gsr1[indr_s])
           IF N_ELEMENTS(indr_s) gt 1 THEN  hbinsrssdv[i] = STDEV([gsr1[indr_s]]) ELSE hbinsrssdv[i] = 0
            hbinsvs[i] = MEDIAN(gsvr1[indr_s])
           IF N_ELEMENTS(indr_s) gt 1 THEN  hbinsrvsdv[i] = STDEV([gsvr1[indr_s]]) ELSE hbinsrvsdv[i] = 0
        ENDIF ELSE BEGIN
            hbinsls[i] = 0
            hbinslssdv[i] = 0
            hbinsrs[i] = 0
            hbinsrssdv[i] = 0
        ENDELSE
        ENDIF
;        oplot,gmatch1[indr].x,gmatch1[indr].z,psym = 3,color = colorsr[i]
;        oplot,ggmatch0[indr].x,ggmatch0[indr].y,psym = 3,color = colorsr[i]
;        stop
    ENDFOR
    oplot,[0,0],[-20,20]
    oplot,[-20,20],[0,0]
    stop
    ind_hot = where((gmatch1.tempg ge maxt) and jvecgg0[*,2] ge 0 and jvecg1[*,2] ge 0); AND coolontime gt 0)

;    plot,gr1,jvecg1[*,2] - jvecgg0[*,2],psym = 3,xtitle = 'Final Radius [kpc]',ytitle = textoidl('L_1 - L_0'),xrange = [0,maxr]
;    oplot,gsr1,jvecs1[*,2] - jvecgs0[*,2],psym = 3,color = 100
;    oplot,[0,maxr],[0,0],color = 50

    plot,ggr0,(jvecg1[*,2] - jvecgg0[*,2])/jvecgg0[*,2],psym = 3,xtitle = 'Initial Radius [kpc]',ytitle = textoidl('L_1 - L_0'),ystyle = 1,xrange = [0,30],/ylog,yrange = [1e-2,1e3];yrange = [-10,20]
    IF NOT KEYWORD_SET(nostarlog) THEN oplot,gsr0,(jvecs1[*,2] - jvecgs0[*,2])/jvecgs0[*,2],psym = 3,color = 60
;    oplot,gsr0[ind_hot],(jvecg1[ind_hot,2] - jvecgg0[ind_hot,2])/jvecgg0[ind_hot,2],psym = 3,color = 240
    oplot,[0,100],[0,0],color = 0
;    oplot,hbins,hbinsl,psym = symcat(14),color = 140,thick = 2,symsize = 2
;    oploterror,hbins,hbinsl,hbinslsdv,errcolor = 140,errthick = 2,color = 140

;    oplot,hbins,hbinslc,psym = symcat(14),color = 100,thick = 2,symsize = 2
;    oploterror,hbins,hbinslc,hbinslsdvc,errcolor = 100,errthick = 2,color = 100

;    oplot,hbins,hbinslh,psym = symcat(14),color = 240,thick = 2,symsize = 2
;    oploterror,hbins,hbinslh,hbinslhsdv,errcolor = 240,errthick = 2,color = 240
    stop

    plot,ggr0,gr1,psym = 3,xtitle = 'Initial Radius [kpc]',ytitle = textoidl('Final Radius [kpc]'),yrange = [0,40],xrange = [0,8],title = title
    IF NOT KEYWORD_SET(nostarlog) THEN oplot,gsr0,gsr1,psym = 3,color = 60
    oplot,ggr0[ind_hot],gr1[ind_hot],psym = 3  ,color = 254 
    oplot,[0,100],[0,100],color = 0
    oplot,hbins[where(hbinsrh ne 0)],hbinsrh[where(hbinsrh ne 0)],psym = symcat(14),color = 240,thick = 2,symsize = 2
    oploterror,hbins[where(hbinsrh ne 0)],hbinsrh[where(hbinsrh ne 0)],hbinsrhsdv[where(hbinsrh ne 0)],errcolor = 240,errthick = 2,color = 240  
    oplot,hbins,hbinsr,psym = symcat(14),color = 60,thick = 2,symsize = 2
    oploterror,hbins,hbinsr,hbinsrsdv,errcolor = 60,errthick = 2,color = 60
    IF NOT KEYWORD_SET(nostarlog) THEN oplot,hbins[where(hbinsrs ne 0)],hbinsrs[where(hbinsrs ne 0)],psym = symcat(14),color = 40,thick = 2,symsize = 2
    IF NOT KEYWORD_SET(nostarlog) THEN oploterror,hbins[where(hbinsrs ne 0)],hbinsrs[where(hbinsrs ne 0)],hbinsrssdv[where(hbinsrs ne 0)],errcolor = 40,errthick = 2,color = 40
    set_plot,'ps'
    formatplot,/outplot
    device,filename = filename + '_'+ steps[0] + '_redist.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=11
    plot,ggr0,gr1,psym = 3,xtitle = 'Initial Radius [kpc]',ytitle = textoidl('Final Radius [kpc]'),yrange = [0,16],xrange = [0,8],title = title
    IF NOT KEYWORD_SET(nostarlog) THEN oplot,gsr0,gsr1,psym = 3,color = 60
    oplot,ggr0[ind_hot],gr1[ind_hot],psym = 3  ,color = 254  
    oplot,[0,100],[0,100],color = 0,thick = 4
    oplot,hbins[where(hbinsrh ne 0)],hbinsrh[where(hbinsrh ne 0)],psym = symcat(14),color = 240,thick = 2,symsize = 1
    oploterror,hbins[where(hbinsrh ne 0)],hbinsrh[where(hbinsrh ne 0)],hbinsrhsdv[where(hbinsrh ne 0)],errcolor = 240,errthick = 6,color = 240  
    oplot,hbins,hbinsr,psym = symcat(14),color = 60,thick = 6,symsize = 1.5
    oploterror,hbins,hbinsr,hbinsrsdv,errcolor = 60,errthick = 6,color = 60,thick = 6
    IF NOT KEYWORD_SET(nostarlog) THEN oplot,hbins[where(hbinsrs ne 0)],hbinsrs[where(hbinsrs ne 0)],psym = symcat(14),color = 40,thick = 6,symsize = 1.5
    IF NOT KEYWORD_SET(nostarlog) THEN oploterror,hbins[where(hbinsrs ne 0)],hbinsrs[where(hbinsrs ne 0)],hbinsrssdv[where(hbinsrs ne 0)],errcolor = 40,errthick = 6,color = 40
    device,/close
    set_plot,'x'
    formatplot
    stop

    plot,ggr0,gvr1,psym = 3,xrange = [0,8],xtitle = 'Initial Radius [kpc]',ytitle = 'Final Radial Velocity',yrange = [-100,200]
    oplot,ggr0[ind_hot],gvr1[ind_hot],psym = 3  ,color = 240
    oplot,[0,100],[0,0],color = 0
    oplot,hbins[where(hbinsvh ne 0)],hbinsvh[where(hbinsvh ne 0)],psym = symcat(14),color = 240,thick = 3,symsize = 2
    oploterror,hbins[where(hbinsvh ne 0)],hbinsvh[where(hbinsvh ne 0)],hbinsvhsdv[where(hbinsvh ne 0)],errcolor = 240,errthick = 3,color = 240,thick = 3
    oplot,hbins[where(hbinsv ne 0)],hbinsv[where(hbinsv ne 0)],psym = symcat(14),color = 60,thick = 1,symsize = 2
    oploterror,hbins[where(hbinsv ne 0)],hbinsv[where(hbinsv ne 0)],hbinsvsdv[where(hbinsv ne 0)],errcolor = 60,errthick = 3,color = 60,thick = 3
    set_plot,'ps'
    formatplot,/outplot
    device,filename = filename + '_'+ steps[0]  + '_redistv.eps',/color,bits_per_pixel= 8,/times,ysize=18,xsize=18
    plot,ggr0,gvr1,psym = 3,xrange = [0,8],xtitle = 'Initial Radius [kpc]',ytitle = 'Final Radial Velocity',yrange = [-100,200]
    oplot,ggr0[ind_hot],gvr1[ind_hot],psym = 3  ,color = 240
    oplot,[0,100],[0,0],color = 0,thick = 6
    oplot,hbins[where(hbinsvh ne 0)],hbinsvh[where(hbinsvh ne 0)],psym = symcat(14),color = 240,thick = 6,symsize = 1.5
    oploterror,hbins[where(hbinsvh ne 0)],hbinsvh[where(hbinsvh ne 0)],hbinsvhsdv[where(hbinsvh ne 0)],errcolor = 240,errthick = 6,color = 240,thick =  6 
    oplot,hbins[where(hbinsv ne 0)],hbinsv[where(hbinsv ne 0)],psym = symcat(14),color = 60,thick = 6,symsize = 1.5
    oploterror,hbins[where(hbinsv ne 0)],hbinsv[where(hbinsv ne 0)],hbinsvsdv[where(hbinsv ne 0)],errcolor = 60,errthick = 6,color = 60,thick = 6
    device,/close
    set_plot,'x'
    formatplot
    stop

    plot,ggr0,gr1 - ggr0,psym = 3,xtitle = 'Initial Radius [kpc]',ytitle = textoidl('Final Radius  - Initial Radius [kpc]'),yrange = [-10,25],xrange = [0,25]
    IF NOT KEYWORD_SET(nostarlog) THEN oplot,gsr0,gsr1 - gsr0,psym = 3,color = 60
 ;   oplot,ggr0[ind_hot],gr1[ind_hot] - ggr0[ind_hot],psym = 3  ,color = 240  
    oplot,[0,100],[0,0],color = 0
    stop

    plot,gmatch1.x,gmatch1.z,psym = 3

    plot,ggr0,gmatch1.tempg - ggmatch0.tempg,psym = 3,xtitle = 'Radius [kpc]',ytitle = textoidl('T_0 - T_1'),/ylog
    stop

    plot,hbins,hbinst,psym = 10,xtitle = 'Initial Radius [Kpc]',ytitle = 'Mass Fraction of Gas Heated by Feedback'  
    stop
ENDIF
  

;    match,iord_gas,star_info.iordergas,temp,ind1_sl;Find the indicies of the starlog file data that represent stars from from gas with the same iord values of those in the file
;    star_info_match = star_info[ind1_sl]
;    match,star_info_match.IORDERSTAR,iord1_star,temp,ind1_star
;    smatch1 = s1[ind1_star]
;    plot,g1[ind1_gas].x,g1[ind1_gas].z,psym = 3,xrange = [-0.003,0.003],yrange = [-0.003,0.003]
;    oplot,s1[ind1_star].x,s1[ind1_star].z,psym = 3,color = 240
;    oplot,g0[ind0_gas].x,g0[ind0_gas].z,psym = 3,color = 40
;    smatch1 = s1[ind1_star]

IF 0 THEN BEGIN
    FOR i = 0, maxr/dr DO BEGIN
        da = 20                 ;0.05
        angmom, gmatch0, jvecg0, lvec0, jmagg0, ltot0
        ind = where(jvecg0[*,2] ge 0.)
;    ind = where(jvecg0[*,2] ge -1e10)
        gmatch0 = gmatch0[ind]
        hg0 = histogram(jvecg0[ind,2], binsize=da, reverse_indices=ri) ;/mean(jvecg0[ind,2])
;hg = histogram(jmagg/mean(jmagg), binsize=da, reverse_indices=ri)
        gmass0 = fltarr(n_elements(hg0))
        for j=0L,n_elements(hg0)-1 do if ri[j+1] gt ri[j] then gmass0[j] = total(gmatch0[ri[ri[j]:ri[j+1]-1]].mass)
;    stop
        
        angmom, gmatch1, jvecg1, lvec1, jmagg1, ltot1
        ind = where(jvecg1[*,2] ge 0.)
;    ind = where(jvecg1[*,2] ge -1e10)
        gmatch1 = gmatch1[ind]
        hg1 = histogram(jvecg1[ind,2], binsize=da, reverse_indices=ri)
;hg = histogram(jmagg/mean(jmagg), binsize=da, reverse_indices=ri)
        gmass1 = fltarr(n_elements(hg1))
        for j=0L,n_elements(hg1)-1 do if ri[j+1] gt ri[j] then gmass1[j] = total(gmatch0[ri[ri[j]:ri[j+1]-1]].mass)
        
        angmom, smatch1, jvecs1, lvec1, jmags1, ltot1
        ind = where(jvecs1[*,2] ge 0.)
;    ind = where(jvecs1[*,2] ge -1e10)
        smatch1 = smatch1[ind]
        hs1 = histogram(jvecs1[ind,2], binsize=da, reverse_indices=ri)
        smass1 = fltarr(n_elements(hs1))
        for j=0L,n_elements(hs1)-1 do if ri[j+1] gt ri[j] then smass1[j] = total(smatch1[ri[ri[j]:ri[j+1]-1]].mass)
        
        
        IF 1 THEN BEGIN
            if(i eq 0) THEN plot, indgen(n_elements(hg0))*da, gmass0/TOTAL(gmass0*da), yrange=[0,0.01], xtickname=xtickname, ytickname=ytickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)',/nodata,title = label, xrange=[0,1200] ;, xrange=[0,3.5], yrange=[0,0.2]
            oplot,indgen(n_elements(hg0))*da, gmass0/TOTAL(gmass0*da),color = color[i]
            if n_elements(hs1) lt n_elements(hg1) then b = n_elements(hg1) else b = n_elements(hs1)
;    oplot, indgen(b)*da, (gmass1/TOTAL(gmass1)+smass1/TOTAL(smass1)), linestyle=2, color = color[i]
            oplot, indgen(n_elements(hg1))*da, (gmass1/TOTAL(gmass1*da)), linestyle=2, color = color[i]
            oplot, indgen(n_elements(hs1))*da, (smass1/TOTAL(smass1*da)), linestyle=3, color = color[i]
        ENDIF
        
        stop
    ENDFOR
    ENDIF
END
