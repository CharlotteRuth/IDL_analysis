;How does the mass of the accreted gas, total gas, and outflowing gas
;from Alyson's Program compare to that found through Sarah L's database

;run from steps directory
;gasmasstrace,'h986.cosmo50cmb.3072g14HBWK','00512','1'

pro gasmasstrace,base,step,halo
infile = base + '.' + step + '.dir/' + base + '.' + step
statfile = infile+'.amiga.stat'
grpfile = infile+'.amiga.grp'
;gtpfile = infile+'.amiga.gtp'
stat = read_stat_struc_amiga(statfile)
grp = read_ascii_array(grpfile)
units = tipsyunits('../' + base + '.param')
rtipsy, infile + '.halo.1', h,g,d,s
znow = (1./h.time)-1.
;rads = sqrt(s.x*s.x+s.y*s.y+s.z*s.z)
;radg = sqrt(g.x*g.x+g.y*g.y+g.z*g.z)
;radd = sqrt(d.x*d.x+d.y*d.y+d.z*d.z)
;cutout_rad = MAX([rads,radg,radd])*units.lengthunit

IF 0 THEN BEGIN
;readcol,base + '.00512.dir/h516_H2_distinct_iord_NOSUB.ascii',sliord_nosub
;readcol,base + '.00512.dir/h516_H2_distinct_iord.ascii',sliord
readcol,base + '.00512.dir/h516_noH2_distinct_iord_NOSUB.ascii',sliord_nosub
readcol,base + '.00512.dir/h516_noH2_distinct_iord.ascii',sliord
ENDIF

IF 1 THEN BEGIN
   inflowz_all  = mrdfits('../grp' + halo + '.accrz.fits')
   inflowm_all  = mrdfits('../grp' + halo + '.mass_at_accr.fits')
   outflowz_all = mrdfits('../grp' + halo + '.outflow_z.fits')
   outflowm_all = mrdfits('../grp' + halo + '.mass_at_outflow.fits')

   inflowz = inflowz_all[where(inflowz_all gt znow)]
   inflowm = inflowm_all[where(inflowz_all gt znow)]
   outflowz = outflowz_all[where(outflowz_all ge znow AND outflowz_all ne 99.0)]
   outflowm = outflowm_all[where(outflowz_all ge znow AND outflowz_all ne 99.0)]
;stop
ENDIF

data = mrdfits(infile + '.grp' + halo + '.allgas.history.fits',1)
data1 = data
data1.x = data.x - (stat[0].xc*1000 - 25000/2)
data1.y = data.y - (stat[0].yc*1000 - 25000/2)
data1.z = data.z - (stat[0].zc*1000 - 25000/2)
radius = sqrt(data1.x*data1.x + data1.y*data1.y + data1.z*data1.z)

inhalo = where(radius le stat[0].rvir, complement = outhalo)

print,'Final Baryonic Mass:   ',total([g.mass,s.mass])*units.massunit
print,'Initial gas mass:      ',3300.0*N_ELEMENTS(data1)
print,'Gas mass outside halo: ',TOTAL(data1[outhalo].mass)
print,'Gas mass inside halo:  ',TOTAL(data1[inhalo].mass)
print,'Percent mass lost:     ',TOTAL(data1[outhalo].mass)/(3300.0*N_ELEMENTS(data1))
stop
end
