;This compares the stellar mass function for the X2X2 and X5X5 runs
;mass_func,'1'

PRO mass_func, outputfile
;msol = 2.362e5
msol = 2.0e17 ; mass of Sum in grams /1e16
dMsolUnit = 1.69875e16 ; Solar mass in system units
kpc = 3.085 ; km per kpc /1e16
grav = 6.67e-23
dKpcUnit = 50000.0 ;Kpc in system units
SYS_AMUCC = 0.0096302124  ;  conversion between system units and amu/c
hubble = 0.70
loadct,39

filename1 = 'Data_files/X2X2randv.dat'
filename2 = 'Data_files/X5X5randv.dat'
;filename3 = 'Data_files/X3X3randv.dat'
;filename4 = 'Data_files/X2X2nfrandv.dat'

readcol,filename1,indicies1,r2001, v2001, s2001, b2001, r11, v11, s11, b11, r21, v21, s21, b21,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
readcol,filename2,indicies2,r2002, v2002, s2002, b2002, r12, v12, s12, b12, r22, v22, s22, b22,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
;readcol,filename3,indicies3,r2003, v2003, s2003, bm2003, r13, v13, s13, bm13, r23, v23, s23, bm23,format = 'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
;readcol,filename4,indicies4,r2004, v2004, s2004, bm2004, r14, v14,
;s14, bm14, r24, v24, s24, bm24,format =
;'A,F,F,F,F,F,F,F,F,F,F,F,F',/silent
m11 = v11^2 * r11 * kpc/grav/msol
m12 = v12^2 * r12 * kpc/grav/msol

zsind1 = WHERE(s11 eq 0)
zsind2 = WHERE(s12 eq 0)
if (zsind1[0] ne -1)then s11[zsind1] = 1;
if (zsind2[0] ne -1)then s12[zsind2] = 1;

zbind1 = WHERE(b11 eq 0)
zbind2 = WHERE(b12 eq 0)
if (zbind1[0] ne -1)then b11[zbind1] = 1;
if (zbind2[0] ne -1)then b12[zbind2] = 1;

log_s1 = ALOG10(s11)
log_s2 = ALOG10(s12)
log_b1 = ALOG10(b11)
log_b2 = ALOG10(b12)
log_m1 = ALOG10(m11)
log_m2 = ALOG10(m12)

mins = 4
maxs = 12
bins = 0.5
s_xaxes = findgen((maxs - mins)/bins)*bins + mins

minb = 4
maxb = 12
binb = 0.5
b_xaxes = findgen((maxb - minb)/binb)*binb + minb

minm = 4
maxm = 12
binm = 0.5
m_xaxes = findgen((maxm - minm)/binm)*binm + minm

s_hist1 = HISTOGRAM(log_s1,binsize = bins, max = maxs, min = mins)
s_hist2 = HISTOGRAM(log_s2,binsize = bins, max = maxs, min = mins)
b_hist1 = HISTOGRAM(log_b1,binsize = binb, max = maxb, min = minb)
b_hist2 = HISTOGRAM(log_b2,binsize = binb, max = maxb, min = minb)
m_hist1 = HISTOGRAM(log_m1,binsize = binm, max = maxm, min = minm)
m_hist2 = HISTOGRAM(log_m2,binsize = binm, max = maxm, min = minm)

set_plot,'ps'
device,filename='s_mfunc.eps',/color,bits_per_pixel=8
plot,s_xaxes,s_hist1,psym = 10,xtitle = "Log Mass (Solar Mass)", ytitle = "Number",title = "Stellar Mass Function"
oplot,s_xaxes,s_hist1,psym = 10, color = 50
oplot,s_xaxes,s_hist2,psym = 10, color = 240
legend,['X2X2','X5X5'],color=[50,240],linestyle=[0,0]
device,/close
;stop

device,filename='b_mfunc.eps',/color,bits_per_pixel=8
plot,b_xaxes,b_hist1,psym = 10,xtitle = "Log Mass (Solar Mass)", ytitle = "Number",title = "Baryonic Mass Function" 
oplot,b_xaxes,b_hist1,psym = 10, color = 50
oplot,b_xaxes,b_hist2,psym = 10, color = 240
legend,['X2X2','X5X5'],color=[50,240],linestyle=[0,0]
device,/close

device,filename='t_mfunc.eps',/color,bits_per_pixel=8
plot,m_xaxes,m_hist1,psym = 10,xtitle = "Log Mass (Solar Mass)", ytitle = "Number",title = "Total Mass Function" 
oplot,m_xaxes,m_hist1,psym = 10, color = 50
oplot,m_xaxes,m_hist2,psym = 10, color = 240
legend,['X2X2','X5X5'],color=[50,240],linestyle=[0,0]
device,/close

device,filename='sm_bm.eps',/color,bits_per_pixel=8
plot,[4,12],[4,12],xtitle = 'Log Baryonic Mass (Solar Mass)',ytitle = 'Log Stellar Mass (Solar Mass)',xrange = [4,12],yrange = [4,12]
oplot,log_b1,log_s1,psym = 2,color = 50
oplot,log_b2,log_s2,psym = 2, color = 240
legend,['X2X2','X5X5'],color=[50,240],linestyle=[0,0],psym = [2,2]
device,/close

device,filename='sm_tm.eps',/color,bits_per_pixel=8
plot,[4,12],[4,12],xtitle = 'Log Total Mass (Solar Mass)',ytitle = 'Log Stellar Mass (Solar Mass)',xrange = [4,12],yrange = [4,12]
oplot,log_m1,log_s1,psym = 2,color = 50
oplot,log_m2,log_s2,psym = 2, color = 240
legend,['X2X2','X5X5'],color=[50,240],linestyle=[0,0],psym = [2,2]
device,/close

set_plot,'x'
;device,filename='bm_tm.eps',/color,bits_per_pixel=8
plot,[4,12],[4,12],xtitle = 'Log Total Mass (Solar Mass)',ytitle = 'Log Baryonic Mass (Solar Mass)',xrange = [4,12],yrange = [4,12],title = 'Halos with no stars'
oplot,log_m1[zsind1],log_b1[zsind1],psym = 2,color = 50
oplot,log_m2[zsind2],log_b2[zsind2],psym = 2, color = 240
legend,['X2X2','X5X5'],color=[50,240],linestyle=[0,0],psym = [2,2]
;device,/close
stop
END
