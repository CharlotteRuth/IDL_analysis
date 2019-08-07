;Very similar to calc_sn but slightly cleaner

 ; .r calc_sn





PRO test_metal_production
dirbase = '/nobackupp8/crchrist/MolecH/' ;PFE
;dirbase = '/home/christensen/Storage2/UW/MolecH/Cosmo/';ymir
dirbase = '/home/christenc/Data/Sims/' ;quirm

steps = ['00024','00048','00072','00096','00120','00144','00168','00192','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']
steps = ['00024','00096','00168','00240','00312','00360','00480','00512']
steps = ['00512']
;steps = ['004096']

;fileroot_short = 'h516.cosmo25cmb.1536g'
;fileroot = 'h516.cosmo25cmb.1536g14HBWK' 
;fileroot_short = 'h516.cosmo25cmb.3072g'
;fileroot = 'h516.cosmo25cmb.3072g14HBWK'
;fileroot_short = 'h603.cosmo50cmb.3072g'
;fileroot = 'h603.cosmo50cmb.3072g14HBWK'
fileroot_short = 'h799.cosmo25cmb.3072g/'
fileroot = 'h799.cosmo25cmb.3072g14HBWK'
;fileroot_short = 'cptmarvel.cosmo25cmb/'
;fileroot = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'

dir = dirbase + fileroot_short + '/' + fileroot
filebase = fileroot + '.'+ steps + '/' + fileroot + '.' + steps
pfile = fileroot + '.param'
;dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/'
;filebase = 'h516.cosmo25cmb.3072g14HBWK.' + steps + '.dir/h516.cosmo25cmb.3072g14HBWK'
;pfile = '../h516.cosmo25cmb.3072g14HBWK.param'

dDelta = 0.000690912 ;g14            ;4.075253e-04
;dDelta = 8.63942e-05 ;Marvel
dDeltaStarForm = 1e6

IF 0 THEN BEGIN
   dir = '/home/christensen/Code/gasoline-radiation/test/onestar/'
   pfile = 'onestar.param'
;steps = ['00020','00031','00032','00033','00034','00035','00036','00037','01000']
;steps = ['00001','00002','00003','00004','00005','00006','00007','00008','00009','00010','00011','00012','00013','00014','00015','00016','00017','00018','00019','00020','00021','00022','00023','00024','00025','00026','00027','00028','00029','00030','00031','00032','00033','00034','00035','00036','00037','00038','00039','00040']
   steps = ['01000','02000','03000','04000','05000','06000','07000','08000','09000','10000']
;steps = ['00100','00200','00300','00400','00500','00600','00700','00800','00900','01000']
;steps = ['00010','00020','00030','00040','00050','00060','00070','00080','00090','00100']
   filebase = 'onestar.' + steps
;filebase = 'onestar.z0.01.' + steps
   dDelta = 7.812500e-04
ENDIF

cd,dir
units = tipsyunits(pfile)
dDelta = dDelta*units.timeunit
IF NOT keyword_set(dDeltaStarForm) THEN dDeltaStarForm = dDelta
totalOxMass = fltarr(n_elements(steps))
totalFeMass = fltarr(n_elements(steps))
totalOxOut = fltarr(n_elements(steps))
totalFeOut = fltarr(n_elements(steps))

FOR i=0,n_elements(steps)- 1 DO BEGIN
    step = steps[i]
    filename = filebase[i]
    rtipsy,filename,h,g,d,s
    readarr,filename + '.FeMassFrac',h,fe_gas,/ascii,part = 'gas',type = 'float'
    readarr,filename + '.FeMassFrac',h,fe_star,/ascii,part = 'star',type = 'float'
    readarr,filename + '.OxMassFrac',h,ox_gas,/ascii,part = 'gas',type = 'float'
    readarr,filename + '.OxMassFrac',h,ox_star,/ascii,part = 'star',type = 'float'
    readarr,filename + '.massform',h,starmass,/ascii,part = 'star',type = 'float'

;sl = rstarlog('../../h516.cosmo25cmb.1536g14HBWK.starlog')
;s.mass = sl.massform ;make sure these line up

    s.tform = s.tform*units.timeunit
    sfmass = starmass*units.massunit
    s.mass = s.mass*units.massunit

    time = max(s.tform)
;time = 0.354154*units.timeunit    
    time = dDelta*float(step)
;stop
    
    which_imf = 0               ;kroupa93, g14
    which_imf = 1  ;kroupa01, marvel
    tmax = time + dDelta
    tmin = time                 ;- dDelta;0 ;time ;time
    tmin = 0
    
    print,''
    print,i,tmin,tmax,time/units.timeunit
    z_outflow_snii = zsnovaii(tmin, tmax, s, sfmass, which_imf)
    z_outflow_snia = zsnovaia(tmin, tmax, s, sfmass)
 
IF pfile eq 'onestar.param' THEN BEGIN
    zstart_ox =  total(sfmass*ox_star)
    zstart_fe =  total(sfmass*fe_star)
ENDIF ELSE BEGIN
    zstart_ox = 0
    zstart_fe = 0
ENDELSE 

    print,''
    print,'Totals'
    print,'Total Ox Mass: ',(total(g.mass*ox_gas)*units.massunit + total(s.mass*ox_star)),', Stars: ',total(s.mass*ox_star),', Gas: ',total(g.mass*ox_gas)*units.massunit
    totalOxMass[i] = (total(g.mass*ox_gas)*units.massunit + total(s.mass*ox_star))
    print,'Total Ox Out: ',z_outflow_snii.oxmassloss + z_outflow_snia.oxmassloss +  zstart_ox,', from SNII: ',z_outflow_snii.oxmassloss,', from SNIa: ',z_outflow_snia.oxmassloss
    totalOxOut[i] = z_outflow_snii.oxmassloss + z_outflow_snia.oxmassloss +  zstart_ox
    print,'Total Fe Mass: ',(total(g.mass*fe_gas)*units.massunit + total(s.mass*fe_star)),', Stars: ',total(s.mass*fe_star),', Gas: ',total(g.mass*fe_gas)*units.massunit
    totalFeMass[i] = (total(g.mass*fe_gas)*units.massunit + total(s.mass*fe_star))
    print,'Total Fe Out: ',z_outflow_snii.femassloss + z_outflow_snia.femassloss + zstart_fe,', from SNII: ',z_outflow_snii.femassloss,', from SNIa: ',z_outflow_snia.femassloss
    totalFeOut[i] = z_outflow_snii.femassloss + z_outflow_snia.femassloss + zstart_fe
ENDFOR
totalzMass =  1.06*totalFeMass + 2.09*totalOxMass
totalzOut = 1.06*totalFeOut + 2.09*totalOxOut

loadct,39
window,0
plot,float(steps)*dDelta,totalzMass
oplot,float(steps)*dDelta,totalzOut,linestyle = 2
oplot,float(steps)*dDelta,totalOxMass,color = 60
oplot,float(steps)*dDelta,totalOxOut,color = 60,linestyle = 2
oplot,float(steps)*dDelta,totalFeMass,color = 254
oplot,float(steps)*dDelta,totalFeOut,color = 254,linestyle = 2

window,1
plot,float(steps)*dDelta,totalzOut/totalzMass
oplot,float(steps)*dDelta,totalOxOut/totalOxMass,color = 60
oplot,float(steps)*dDelta,totalFeOut/totalFeMass,color = 254
stop
END

;print,(total(g.mass*femassfrac)*units.massunit + total(s.mass*femassfrac))
;      280469.

;print,(total(g.mass*oxmassfrac)*units.massunit + total(s.mass*oxmassfrac))
;  2.37371e+06
