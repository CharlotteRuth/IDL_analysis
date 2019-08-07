;Charlotte Christensen
;8/22/18
;Written to examine what the properties of star forming gas are

units = tipsyunits('storm.cosmo25cmb.4096g5HbwK1BH.param')
;/home/christenc/Code/IDL/IDL_analysis/idl4tipsy/rtipsy.pro
rtipsy,'storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096',h,g,d,s
sl = rstarlog('storm.cosmo25cmb.4096g5HbwK1BH.starlog.temp',/molecularH,/big)
read_tipsy_arr,'storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096.iord',h,siord,part = 'star',type = 'long'

match,siord,sl.iorderstar,ind1,ind2
siord = siord[ind1]
sl = sl[ind2]
