filename = 'h603.cosmo50cmb.3072g14HBWK.00268'
filename1 = '../h603.cosmo50cmb.3072g14HBWK.00272.dir/h603.cosmo50cmb.3072g14HBWK.00272'
filename0 = '../h603.cosmo50cmb.3072g14HBWK.00264.dir/h603.cosmo50cmb.3072g14HBWK.00264'

filename = 'h603.cosmo50cmb.3072g14HBWK.00476'
filename1 = '../h603.cosmo50cmb.3072g14HBWK.00480.dir/h603.cosmo50cmb.3072g14HBWK.00480'
filename0 = '../h603.cosmo50cmb.3072g14HBWK.00472.dir/h603.cosmo50cmb.3072g14HBWK.00472'

filenamelast = '../h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512'
slfile = '../../h603.cosmo50cmb.3072g14HBWK.starlog.fits'
sl = mrdfits(slfile,1)

rtipsy,filename,h,g,d,s,/justhead
rtipsy,filename0,h0,g0,d0,s0,/justhead
read_tipsy_arr,filename0 + '.iord',h0,iord0,type = 'long'
iord0_g = iord0[0:h0.ngas - 1]
iord0_d = iord0[h0.ngas:h0.ngas + h0.ndark -  1]
iord0_s = iord0[h0.ngas + h0.ndark:h0.ngas + h0.ndark + h0.nstar -  1]

rtipsy,filename1,h1,g1,d1,s1,/justhead
read_tipsy_arr,filename1 + '.iord',h1,iord1,type = 'long'
iord1_g = iord1[0:h1.ngas - 1]
iord1_d = iord1[h1.ngas:h1.ngas + h1.ndark -  1]
iord1_s = iord1[h1.ngas + h1.ndark:h1.ngas + h1.ndark +h1.nstar -  1]

iord_d = iord1_d
iord_s = iord1_s[0:h.nstar - 1]

match2,iord0_g,iord1_g,ind0,ind1
iordg0_missing = iord0_g[where(ind0 EQ -1)] ;Gas particles that dissapear between first and last step
iord_sformed1 = iord1_s[h.nstar:h1.nstar - 1] ;Stars formed bewteen main step and last step
match2,sl.iorderstar,iord_sformed1,ind2,ind3
sl_formed = sl[where(ind2 NE -1)]
match2,sl_formed.iordergas,iordg0_missing,ind4,ind5 ;"Missing" gas particles stars are formed from between main step and last step
iordg_left = iordg0_missing[where(ind5 NE -1)]
iordg_left = iordg_left[uniq(iordg_left,sort(iordg_left))]
match2,[iordg_left,iord1_g],iord0_g,ind6,ind7 ;Resort gas particles into the correct order
;iordg_left = iordg0_missing[where(ind5 EQ -1)]
;iordg_left = iordg_left[uniq(iordg_left,sort(iordg_left))]
;match2,[iordg_left,iord1_g],iord0_g,ind6,ind7
iord_g = iord1_g[where(ind7 NE -1)]
iord = [iord_g,iord_d,iord_s]

openw,lunhi,filename + '.iord',/get_lun
printf,lunhi,h.n
for j=0L,h.n-1 do printf,lunhi,iord[j]
close,lunhi

rtipsy,filenamelast,hlast,glast,dlast,slast,/justhead
read_tipsy_arr,filenamelast + '.iord',hlast,iordlast,type = 'long',part = 'star'
read_tipsy_arr,filenamelast + '.igasorder',hlast,igasorderlast,type = 'long',part = 'star'
read_tipsy_arr,filenamelast + '.massform',hlast,mflast,type = 'float',part = 'star'
read_tipsy_arr,filenamelast + '.timeform',hlast,tflast,type = 'float',part = 'star'

match2,iord,iordlast,ind0,ind1
igasorder = igasorderlast[where(ind0 NE -1)]
mf = mflast[where(ind0 NE -1)]
tf = tflast[where(ind0 NE -1)]
igasorder = [fltarr(h.ngas),fltarr(h.ndark),igasorder]
mf = [fltarr(h.ngas),fltarr(h.ndark),mf]
tf = [fltarr(h.ngas),fltarr(h.ndark),tf]

openw,lunhi,filename + '.iordergas',/get_lun
printf,lunhi,h.n
for j=0L,h.n-1 do printf,lunhi,igasorder[j]
close,lunhi

openw,lunhi,filename + '.massform',/get_lun
printf,lunhi,h.n
for j=0L,h.n-1 do printf,lunhi,mf[j]
close,lunhi

openw,lunhi,filename + '.timeform',/get_lun
printf,lunhi,h.n
for j=0L,h.n-1 do printf,lunhi,tf[j]
close,lunhi
