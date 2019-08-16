;Very similar to calc_sn but slightly cleaner

 ; .r calc_sn

FUNCTION int_kroupa01,mass1,mass2,nstar = nstar,mod_exp = mod_exp
;Integratees the IMF to give total number of stars between bounding
;masses, mass1 and mass2, where mass1 is the greater mass

IF n_elements(mass1) NE n_elements(mass2) THEN BEGIN
    print,'mass1 and mass2 must have the same number of elements'
    RETURN,0*mass1
ENDIF

bmasses = [0.08,0.5,100] ;bounding masses and mass of break in IMF
exp = [-0.3,-1.3] ;exponents of IMF to produce mass of stars in given range. IMF to produce # will have exp-1
a = [2.0*alog(10.0),alog(10.0)]*0.22038 ;/14.4981 ;coefficients for IMF
IF keyword_set(nstar) THEN exp = exp-1  ;so that it is no longer in terms of mass but is in number
IF keyword_set(mod_exp) THEN exp = exp + mod_exp ;Enables one to include a function of stellar mass (for instance, oxygen mass released) inside of the integral

IF (where(mass1 GT bmasses[2]))[0] NE -1 THEN mass1[where(mass1 GT bmasses[2])] = bmasses[2]
IF (where(mass2 GT bmasses[2]))[0] NE -1 THEN mass2[where(mass2 GT bmasses[2])] = bmasses[2]
IF (where(mass1 LT bmasses[0]))[0] NE -1 THEN mass1[where(mass1 LT bmasses[0])] = bmasses[0]
IF (where(mass2 LT bmasses[0]))[0] NE -1 THEN mass2[where(mass2 LT bmasses[0])] = bmasses[0]

cumM1 = fltarr(n_elements(mass1)) ;Create an array to hold the total mass between mass1 and the maximum mass

ind1 = where(mass1 GT bmasses[1] AND mass1 LE bmasses[2]) ;Anything that is in the high-mass part of the IMF
ind2 = where(mass1 GE bmasses[0] AND mass1 LE bmasses[1]) ;Anything that is in the low-mass part of the IMF

IF (ind1[0] NE -1) THEN cumM1[ind1] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - mass1[ind1]^(exp[1] + 1))
IF (ind2[0] NE -1) THEN cumM1[ind2] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - bmasses[1]^(exp[1] + 1)) + $
                                      a[0]/(exp[0] + 1)*(bmasses[1]^(exp[0] + 0) - mass1[ind2]^(exp[0] + 1))

cumM2 = fltarr(n_elements(mass1)) ;Create an array to hold the total mass between mass2 and the maximum mass

ind1 = where(mass2 GT bmasses[1] AND mass2 LE bmasses[2]) ;Anything that is in the high-mass part of the IMF
ind2 = where(mass2 GE bmasses[0] AND mass2 LE bmasses[1]) ;Anything that is in the low-mass part of the IMF

IF (ind1[0] NE -1) THEN cumM2[ind1] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - mass2[ind1]^(exp[1] + 1))
IF (ind2[0] NE -1) THEN cumM2[ind2] = a[1]/(exp[1] + 1)*(bmasses[2]^(exp[1] + 1) - bmasses[1]^(exp[1] + 1)) + $
                                      a[0]/(exp[0] + 1)*(bmasses[1]^(exp[0] + 0) - mass2[ind2]^(exp[0] + 1))

RETURN,cumM2 - cumM1
END




PRO test_metal_production
dirbase = '/nobackupp8/crchrist/MolecH/' ;PFE
;dirbase = '/home/christensen/Storage2/UW/MolecH/Cosmo/';ymir
dirbase = '/home/christenc/Data/Sims/' ;quirm

steps = ['00024','00048','00072','00096','00120','00144','00168','00192','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']
steps = ['00024','00096','00168','00240','00312','00360','00480','00512']
steps = ['00512']
;fileroot_short = 'h516.cosmo25cmb.1536g'
;fileroot = 'h516.cosmo25cmb.1536g14HBWK' 
;fileroot_short = 'h516.cosmo25cmb.3072g'
;fileroot = 'h516.cosmo25cmb.3072g14HBWK'
;fileroot_short = 'h603.cosmo50cmb.3072g'
;fileroot = 'h603.cosmo50cmb.3072g14HBWK'
fileroot_short = 'h799.cosmo25cmb.3072g'
fileroot = 'h799.cosmo25cmb.3072g14HBWK'
dDelta = 0.000690912            ;4.075253e-04
which_imf = 0                   ;kroupa93
changa = 0

IF 0 THEN BEGIN
   steps = ['004096']
   fileroot_short = 'cptmarvel.cosmo25cmb'
   fileroot = 'cptmarvel.cosmo25cmb.4096g5HbwK1BH'
   dDelta = 8.63942e-05         ;Marvel
   which_imf = 1                ;kroupa01
   changa = 1
ENDIF
   
dDeltaStarForm = 1e6

dir = dirbase + fileroot_short + '/' + fileroot
filebase = fileroot + '.'+ steps + '/' + fileroot + '.' + steps
pfile = fileroot + '.param'
;dir = '/home/christensen/Storage2/UW/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/'
;filebase = 'h516.cosmo25cmb.3072g14HBWK.' + steps + '.dir/h516.cosmo25cmb.3072g14HBWK'
;pfile = '../h516.cosmo25cmb.3072g14HBWK.param'

IF 0 THEN BEGIN
   dir = '/home/christensen/Code/gasoline-radiation/test/onestar/'
   dir = '/home/christenc/Code/changa_uw/testonestar/'
   pfile = 'onestar.param'
   steps = ['01000','02000','03000','04000','05000','06000','07000','08000','09000','10000']
   steps = ['000130']
   filebase = 'onestar.' + steps
   dDelta = 7.812500e-04
   which_imf = 1
   changa = 1
ENDIF

cd,dir
units = tipsyunits(pfile)
dDelta = dDelta*units.timeunit
print,dDelta*units.timeunit
IF keyword_set(dDeltaStarForm) THEN dDeltaStarForm = dDelta
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
    
    tmax = time + dDelta
    tmin = time                 ;- dDelta;0 ;time ;time
    tmin = 0
    
    print,''
    print,i,tmin,tmax,time/units.timeunit
    z_outflow_snii = zsnovaii(tmin, tmax, s[where(s.tform ge 0)], sfmass[where(s.tform ge 0)], which_imf)
    z_outflow_snia = zsnovaia(tmin, tmax, s[where(s.tform ge 0)], sfmass[where(s.tform ge 0)], changa = changa)
 
IF pfile eq 'onestar.param' THEN BEGIN
    zstart_ox =  total(sfmass*ox_star)
    zstart_fe =  total(sfmass*fe_star)
ENDIF ELSE BEGIN
    zstart_ox = 0
    zstart_fe = 0
ENDELSE 

    print,''
    print,'Totals'
    print,'Total Ox Mass: ',(total(g.mass*ox_gas)*units.massunit + total(s[where(s.tform gt 0)].mass*ox_star[where(s.tform gt 0)])),', Stars: ',total(s[where(s.tform gt 0)].mass*ox_star[where(s.tform gt 0)]),', Gas: ',total(g.mass*ox_gas)*units.massunit
    totalOxMass[i] = (total(g.mass*ox_gas)*units.massunit + total(s[where(s.tform gt 0)].mass*ox_star[where(s.tform gt 0)]))
    print,'Total Ox Out: ',z_outflow_snii.oxmassloss + z_outflow_snia.oxmassloss +  zstart_ox,', from SNII: ',z_outflow_snii.oxmassloss,', from SNIa: ',z_outflow_snia.oxmassloss
    totalOxOut[i] = z_outflow_snii.oxmassloss + z_outflow_snia.oxmassloss +  zstart_ox
    print,'Total Fe Mass: ',(total(g.mass*fe_gas)*units.massunit + total(s[where(s.tform gt 0)].mass*fe_star[where(s.tform gt 0)])),', Stars: ',total(s[where(s.tform gt 0)].mass*fe_star[where(s.tform gt 0)]),', Gas: ',total(g.mass*fe_gas)*units.massunit
    totalFeMass[i] = (total(g.mass*fe_gas)*units.massunit + total(s[where(s.tform gt 0)].mass*fe_star[where(s.tform gt 0)]))
    print,'Total Fe Out: ',z_outflow_snii.femassloss + z_outflow_snia.femassloss + zstart_fe,', from SNII: ',z_outflow_snii.femassloss,', from SNIa: ',z_outflow_snia.femassloss
    totalFeOut[i] = z_outflow_snii.femassloss + z_outflow_snia.femassloss + zstart_fe
ENDFOR
totalzMass =  1.06*totalFeMass + 2.09*totalOxMass
totalzOut = 1.06*totalFeOut + 2.09*totalOxOut

IF n_elements(steps) > 1 THEN BEGIN
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
ENDIF

stop

END

;print,(total(g.mass*femassfrac)*units.massunit + total(s.mass*femassfrac))
;      280469.

;print,(total(g.mass*oxmassfrac)*units.massunit + total(s.mass*oxmassfrac))
;  2.37371e+06
