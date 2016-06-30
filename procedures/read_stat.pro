pro read_stat, FILE = file, NEWCENTERS = newcenters

common halo_stats, group, npart, vc, vc_half, radius, r_vc_max, r_m_half, xcenter, ycenter, zcenter, xVcm, yVcm, zVcm, M_star, m_tot, m_gas
common bound_center, xbound, ybound, zbound

mass_unit = 13.6*10^16.
v_unit = 2418.5
length_unit = 1d5

readcol, file, group, member, m_tot, m_gas, m_star, vc_max, vc_half, vc_outer, r_vc_max, r_m_half, r_outer, v_disp, xcenter, ycenter, zcenter, xVcm, yVcm, zVcm, xbound, ybound, zbound, /silent

npart = member
vc = vc_max

if keyword_set(newcenters) then begin
 tag = 'stat'
 locate = STRPOS(file,tag)
 newtag = 'centers'
 center_file = STRMID(file,0,locate)+newtag
 readcol, center_file, xcenter_new, ycenter_new, zcenter_new
; print, n_elements(xcenter_new), n_elements(xcenter)

 for i = 0, n_elements(xcenter_new)-1 do begin
  if m_star[i] gt 0 then print, (xcenter_new[i] - xcenter[i])*1d5, (ycenter_new[i] - ycenter[i])*1d5, (zcenter_new[i] - zcenter[i])*1d5, r_outer[i]*1d5
 endfor
endif

radius = r_outer

RETURN
END
