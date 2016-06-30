function read_stat_struc_AMIGA, file
;readcol, file, group, member, ngas, nstar, ndark, m_tot, R_vir, $
;  m_gas, m_star, m_dark, vc_max, r_vc_max, v_disp,  xc, yc, zc, $
;  xVcm, yVcm, zVcm, origID, /silent

;if (N_ELEMENTS(group eq 0) gt 1 )then begin
;    readcol, file, group, member, ngas, nstar, ndark, m_tot, R_vir, $
;      m_gas, m_star, m_dark, vc_max, r_vc_max, v_disp,  xc, yc, zc, $
;      xVcm, yVcm, zVcm, origID, satellite, /silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A,A)'
    readcol, file, group, member, ngas, nstar, ndark, m_tot, R_vir, $
      m_gas, m_star, m_dark, vc_max, r_vc_max, v_disp,  xc, yc, zc, $
      xVcm, yVcm, zVcm, contam, satellite, /silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A,A)'
;endif ELSE satellite = STRARR(N_ELEMENTS(group)) 

if (N_ELEMENTS(group) eq 0 OR satellite[0] eq 'false' OR satellite[0] eq 'false?' OR MAX(group) eq 0) then begin
;    readcol, file, group, member, ngas, nstar, ndark, m_tot, R_vir, $
;      m_gas, m_star, m_dark, vc_max, r_vc_max, v_disp,  xc, yc, zc, $
;      xVcm, yVcm, zVcm, origID, /silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A)'
    readcol, file, group, member, ngas, nstar, ndark, m_tot, R_vir, $
      m_gas, m_star, m_dark, vc_max, r_vc_max, v_disp,  xc, yc, zc, $
      xVcm, yVcm, zVcm, contam,/silent, format ='(F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,f,F,F,F,A)'
    satellite = STRARR(N_ELEMENTS(group)) 
endif

;true = where(false ne 'false' AND false ne 'false?')
;group = group[true]
;member = member[true]
;ngas = ngas[true] 
;nstar = nstar[true] 
;ndark = ndark[true]
;m_tot = m_tot[true] 
;R_vir = R_vir[true]
;m_gas = m_gas[true]
;m_star =m_star[true]
;m_dark =m_dark[true]
;vc_max =vc_max[true]
;r_vc_max =r_vc_max[true]
;v_disp =v_disp[true]
;xc =xc[true]
;yc =yc[true]
;zc =zc[true]
;xVcm =xVcm[true]
;yVcm =yVcm[true]
;zVcm =zVcm[true]
;contam =contam[true]
;false =false[true]

;stop

halos = replicate({group:0, npart:0L, m_tot:0., m_gas:0., m_star:0., m_dark:0., $
                  vc:0., r_vc_max:0., ngas:0L, nstar:0L, ndark:0L, rvir:0., $
                  v_disp:0., xc:0., yc:0., zc:0., xVcm:0., $
                   yVcm:0., zVcm:0.,sat:' '}, n_elements(group))

halos.vc = REFORM(vc_max)
halos.group = REFORM(group)
halos.npart = REFORM(member)
halos.ngas = REFORM(ngas)
halos.nstar = REFORM(nstar)
halos.ndark = REFORM(ndark)
halos.m_tot = REFORM(m_tot)
halos.m_gas = REFORM(m_gas)
halos.m_star = REFORM(m_star)
halos.m_dark = REFORM(m_dark)
halos.rvir = REFORM(R_vir)
halos.r_vc_max = REFORM(r_vc_max)
halos.v_disp = REFORM(v_disp)
halos.xc = REFORM(xc)
halos.xVcm = REFORM(xVcm)
halos.yc = REFORM(yc)
halos.yVcm = REFORM(yVcm)
halos.zc = REFORM(zc)
halos.zVcm = REFORM(zVcm)
halos.sat = REFORM(satellite)
RETURN, halos
END
