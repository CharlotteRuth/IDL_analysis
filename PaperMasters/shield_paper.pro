;Run pynbody to find halos
;
;import pynbody
;s = pynbody.load('h516.cosmo25cmb.3072g1HBWK.000068')
;h = s.halos()

dir = '/home/christensen/Storage1/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1HBWKS.00444/'
pfile = 'h516.cosmo25cmb.3072g1HBWKS.00444.param'
tfile = 'h516.cosmo25cmb.3072g1HBWKS.000068'

dir = '/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWK.00444/'
pfile = 'h516.cosmo25cmb.3072g1HBWK.00444.param'
tfile = 'h516.cosmo25cmb.3072g1HBWK.000068'

haloid = '1'
units = tipsyunits(pfile)
ahf_grp_stat,dir+tfile,boxsize = units.LENGTHUNIT*units.h,munit = units.massunit,vunit = units.vunit
outarray = ['HI','H2','iord','FeMassFrac','OxMassFrac','HeI','HeII'] 
tipsysatshi,tfile,1,units.lengthunit,units.massunit,/cutout_rad,outarray = outarray
print,'ln -s '+tfile+'.halo.'+haloid+'.std '+tfile+'.halo.'+haloid
stop
writehimacro,tfile+'.halo.'+haloid,pfile =pfile,radius = 24,/dwarfv,/molecularH
print,'Open tipsy'
print,'run: readmacro '+tfile+'.halo.'+haloid+'HI.macro'
print,'macro makeHIcube'
print,'Open con_script.analyzeHIcubes and make sure it looks like the below'
stop
analyzeHIcubes,dir,tfile+'.halo.'+haloid,angle = 90,/dwarf,/physical_coord
analyzeHIcubes,dir,tfile+'.halo.'+haloid,angle = 45,/dwarf,/physical_coord
schmidtlaw_res_obs,tfile,pfile,extno = 16,/verbose,/useH2,angle = 45
schmidtlaw_res_obs_master_out,['/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWKS.00444/schmidtlaw_res_obs_all_Ha0.74626866.dat','/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWK.00444/schmidtlaw_res_obs_all0.750000.dat'],color = [254,60],/formatthick,key = ['S','H2'],outplot = '~/h516'

;----------------------------------------------
step = '000092'

dirH2 = '/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWK/'
pfileH2 = 'h516.cosmo25cmb.3072g1HBWK.param'
tfileH2 = 'h516.cosmo25cmb.3072g1HBWK.' + step

dirS = '/home/christenc/Storage/Cosmo/h516.cosmo25cmb/h516.cosmo25cmb.3072g1HBWKS/'
pfileS = 'h516.cosmo25cmb.3072g1HBWKS.param' 
tfileS = 'h516.cosmo25cmb.3072g1HBWKS.' + step 
