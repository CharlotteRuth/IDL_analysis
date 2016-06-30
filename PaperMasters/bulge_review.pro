;outplot = '/home/christensen/Plots/fbmass/bulge'
;outplot = '/home/christensen/Plots/fbmass/bulge_thick'
;outplot = '/home6/crchrist/Plots/bulge'
;outplot = '/home6/crchrist/Plots/bulge_thick'

PRO bulge_review,outplot = outplot
formatplot,outplot = outplot,thick = formatthick
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN bgcolor = 255 ELSE bgcolor = 0
loadct,39
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' $
ELSE IF (strcmp(hostname, 'bridge', 6) OR strcmp(hostname, 'pfe', 3)) THEN prefix = '/nobackupp8/crchrist/MolecH/' $
ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

plot_sfh = 0
plot_ml = 0
plot_cm = 0
find_central_mass_history = 0
plot_stellar_dist = 0
follow_central_sf = 1

dir799 = prefix + 'h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/'
file799 = 'h799.cosmo25cmb.3072g14HBWK'
key799 = 'h799'
;116,124,152,172
dir516 = prefix + 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/'
file516 = 'h516.cosmo25cmb.3072g14HBWK'
key516 = 'h516'
dir986 = prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/'
file986 = 'h986.cosmo50cmb.3072g14HBWK'
key986 = 'h986'
dir603 = prefix + 'h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/'
file603 = 'h603.cosmo50cmb.3072g14HBWK'
key603 = 'h603'
dir285 = prefix + 'h285.cosmo50cmb.3072g/h285.cosmo50cmb.3072g14HMbwK/'
file285 = 'h285.cosmo50cmb.3072g14HMbwK'
key285 = 'h285'
dir258 = prefix + 'h258.cosmo50cmb.3072g/h258.cosmo50cmb.3072g14HMbwK/'
file258 = 'h258.cosmo50cmb.3072g14HMbwK'
key258 = 'h258'
dir239 = prefix + 'h239.cosmo50cmb.3072g/h239.cosmo50cmb.3072g14HMbwK/'
file239 = 'h239.cosmo50cmb.3072g14HMbwK'
key239 = 'h239'

dirs   = [dir799, dir516, dir986, dir986, dir986, dir603, dir285, dir258, dir239, dir239]
files  = [file799,file516,file986,file986,file986,file603,file285,file258,file239, file239]
keys   = [key799, key516, key986, key986+'1', key986, key603, key285, key258, key239, key239+'1']
step1  = ['00124','00224','00124','00244','00100','00212','00104','00192','00204','00216']
step2  = ['00172','00268','00168','00284','00140','00260','00168','00252','00260','00288']
step0  = ['00512','00512','00512','00512','00512','00512','00512','00512','00512','00512']
primid = ['1',    '1',    '1',    '1',    '4',    '1',    '1',    '1',    '1',    '4']
secid  = ['4',    '3',    '5',    '4',    '7',    '3',    '2',    '2',    '2',    '14']
finid  = ['1',    '1',    '1',    '1',    '4',    '1',    '1',    '1',    '1',    '2']
z0id   = ['1',    '1',    '1',    '1',    '3',    '1',    '1',    '1',    '1',    '2']
z0mass = textoidl(['2.3 \times 10^{10}','3.8 \times 10^{10}','1.8 \times 10^{11}','1.8 \times 10^{11}','3.8 \times 10^{10}','3.4 \times 10^{11}','8.8 \times 10^{11}','7.7 \times 10^{11}','9.1 \times 10^{11}','4.3 \times 10^{10}']) + 'M' + sunsymbol()
massloc = [0,0,0,0,1,1,1,1,1]
massloc = [1,1,1,1,1,0,0,0,0,0]
outname = ['h799.1','h516.1','h986.1_1','h986.1_2','h986.3','h603.1','h285.1','h258.1','h239.1','h239.2']

IF find_central_mass_history THEN $
   FOR i = 0, n_elements(dirs) - 1 DO central_mass_history,dir = dirs[i]

IF plot_sfh THEN $
   FOR i = 0, n_elements(dirs) - 1 DO bulge_sfh,dirs[i],files[i],haloid = z0id[i],step0 = step0[i],step1 = step1[i],step2 = step2[i],key = keys[i],outplot = outplot,massloc = massloc[i],label = z0mass[i]

IF plot_ml THEN $
   plot_massload,dirs,files,z0id,outplot = outplot

IF plot_cm THEN $
   plot_cmass,dirs,files,z0id,step1,step2,outplot = outplot

IF plot_stellar_dist THEN $
   FOR i = 0, n_elements(dirs) - 1 DO comp_stellar_dist,dirs[i],files[i],z0id[i],step1[i],primid[i],step2[i],finid[i],outplot = outplot
;IF plot_stellar_dist THEN $
;   FOR i = 0, n_elements(dirs) - 1 DO comp_stellar_dist,dirs[i],files[i],step1[i],primid[i],step0[i],z0id[i]


IF follow_central_sf THEN BEGIN
   valid = [0,1,2,4,6,7,8,9]
   FOR i = 1, n_elements(valid) - 1 DO BEGIN
      IF keyword_set(outplot) THEN outplot_step = '~/Plots/' + outname[valid[i]]
      follow_central_sf,files[valid[i]],step1[valid[i]],step2[valid[i]],dir = dirs[valid[i]],finalid = z0id[valid[i]],outplot = outplot_step
   ENDFOR
ENDIF

END
