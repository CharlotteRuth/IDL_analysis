;120 -- 135

PRO exp,x,A,F,pder
;x is the radius, 
;A[0]= effective radius, 
;A[1] = surface brightness at the effective radius
F = A[1]*EXP((-1.0)*X/A[0])
d1 = EXP((-1.0)*X/A[0])
d2 = A[1]*EXP(x/A[0])*x*A[0]^(-2.0)
pder =[[d1],[d2]]
END

function S_H2,yH2,h
omega_H2 = 0.2
x = yH2*h/5d14
return, (1 - omega_H2)/(1 + x)^2 + omega_H2/SQRT(1 + x)*exp(-0.00085*SQRT(1 + x))
end

function S_d,yHI,yH2,h;,z
ZSOLAR = 0.0130215
sigma_d = 2d-21
;return, exp(-1.0*sigma_d*z/ZSOLAR*(yHI*h + 2.0*yH2*h))
return, exp(-1.0*sigma_d*(yHI*h + 2.0*yH2*h))
end

pro rotCurve_12M
prefix = "/astro/net/scratch2/christensen/MolecH/12M/"

base = "/astro/net/scratch2/christensen/MolecH/11M/Disk_Iso_1e5_zsol/uvtest/MW_disk.00009"

keys = ['1E5 Particles, psource']
;legend_t = '1E6 Particles'

msol_per_sysmass = 1.36e17
kpc_per_syslength = 1e5

rotcurve,[base],msol_per_sysmass,kpc_per_syslength,keys =keys,/color,redshift = 0.0001
IF (NOT KEYWORD_SET(last)) THEN last = 300
IF (NOT KEYWORD_SET(dt)) THEN dt = 1
IF (NOT KEYWORD_SET(first)) THEN start = 1 ELSE start = first  ;10/dt + 1

FOR i = start/dt, last/dt DO BEGIN 
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    filenametip = prefix+base+'.'+step+'.'+STRTRIM(cosmo,2)+'.std'
    filename = prefix+base+"."+step
ENDFOR

end

pro rotCurve_h2003
halo = 1
msol_per_sysmass       = 2.310e15
kpc_per_syslength        = 25000.

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
rotcurve,prefix+'h2003.cosmo25cmb.4096g/h2003.cosmo25cmb.4096g14HBK/steps/h2003.cosmo25cmb.4096g14HBK.00512.dir/h2003.cosmo25cmb.4096g14HBK.00512.halo.1',msol_per_sysmass,kpc_per_syslength
end

pro rotCurve_h516
halo = 1
msol_per_sysmass       = 2.310e15
kpc_per_syslength        = 25000.

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

;base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00120.dir/h516.cosmo25cmb.1536g14HBWK.00120.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noJeans/h516.cosmo25cmb.1536g3HBWK_noJeans.00120.dir/h516.cosmo25cmb.1536g3HBWK_noJeans.00120.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00120.dir/h516.cosmo25cmb.1536g1MBWK.00120.halo.1']
;keys = ['new H2','old H2','no H2']
;outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.1536g3H_14H_MBWK.00120'

;base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir/h516.cosmo25cmb.1536g1MBWK.00512.halo.1']
;keys = ['new H2','no H2']
;outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.1536g14H_MBWK.00512'

;base = ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00448.dir/h516.cosmo25cmb.2304g14HBWK.00448.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00432.dir/h516.cosmo25cmb.1536g14HBWK.00432.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noJeans/h516.cosmo25cmb.1536g3HBWK_noJeans.00432.dir/h516.cosmo25cmb.1536g3HBWK_noJeans.00432.halo.1', $
;        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00432.dir/h516.cosmo25cmb.1536g1MBWK.00432.halo.1']
;keys = ['new H2, mr','new H2, lr','old H2','no H2 hr','no H2 lr']
;outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.2304_1536g3H_14H_MBWK.00432'

;base = ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00144.dir/h516.cosmo25cmb.2304g14HBWK.00144.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00144.dir/h516.cosmo25cmb.1536g14HBWK.00144.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g15HBWK/steps/h516.cosmo25cmb.1536g15HBWK.00144.dir/h516.cosmo25cmb.1536g15HBWK.00144.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noJeans/h516.cosmo25cmb.1536g3HBWK_noJeans.00144.dir/h516.cosmo25cmb.1536g3HBWK_noJeans.00144.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00144.dir/h516.cosmo25cmb.1536g1MBWK.00144.halo.1']
;keys = ['new H2, hr','new H2, lr','new H2, high T','old H2','no H2']
;outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.2304_1536g3H_14H_MBWK.00144'

;base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00144.dir/h516.cosmo25cmb.1536g14HBWK.00144.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g15HBWK/steps/h516.cosmo25cmb.1536g15HBWK.00144.dir/h516.cosmo25cmb.1536g15HBWK.00144.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noJeans/h516.cosmo25cmb.1536g3HBWK_noJeans.00144.dir/h516.cosmo25cmb.1536g3HBWK_noJeans.00144.halo.1', $
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00144.dir/h516.cosmo25cmb.1536g1MBWK.00144.halo.1']
;keys = ['new H2, lr','new H2, high T','old H2','no H2']
;outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.2304_1536g3H_14H_MBWK.00144'

base = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/h516.cosmo25cmb.2304g14HBWK.00512.halo.1']
keys = ['no '+textoidl('H_2'), textoidl('H_2')]
outfile = '~/h516.cosmo25cmb.paper'

base = ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00444.dir/h516.cosmo25cmb.2304g14HBWK.00444.halo.1',$
        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00444.dir/h516.cosmo25cmb.3072g14HBWK.00444.halo.1']
keys = ['old '+textoidl('H_2'), 'new '+ textoidl('H_2')]

base = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
keys = ['no '+textoidl('H_2'), textoidl('H_2')]
outfile = '~/h516.cosmo25cmb.paper'

base = [ 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
keys = [textoidl('H_2')]
outfile = '~/h516.cosmo25cmb.paper'

base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.00512.dir/h516.cosmo25cmb.1536g2MBWK.00512.halo.1', $
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g16HBWK/steps/h516.cosmo25cmb.1536g16HBWK.00512.dir/h516.cosmo25cmb.1536g16HBWK.00512.halo.1',$
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.halo.1',$
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g17HBWK/steps/h516.cosmo25cmb.1536g17HBWK.00512.dir/h516.cosmo25cmb.1536g17HBWK.00512.halo.1']
keys = ['no '+textoidl('H_2'), 'no '+textoidl('H_2 SF'), textoidl('H_2'),textoidl('H_2, Gnedin')]
outfile = prefix + 'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536gSF'

base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.00512.dir/h516.cosmo25cmb.1536g2MBWK.00512.halo.1', $
        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/steps/h516.cosmo25cmb.2304g2MBWK.00512.dir/h516.cosmo25cmb.2304g2MBWK.00512.halo.1',$
         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g16HBWK/steps/h516.cosmo25cmb.2304g16HBWK.00512.dir/h516.cosmo25cmb.2304g16HBWK.00512.halo.1',$
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.halo.1',$
        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
keys = ['no '+textoidl('H_2, lr'),'no '+textoidl('H_2, mr'),'no '+textoidl('H_2, hr'), 'no '+textoidl('H_2 SF, mr'), textoidl('H_2, lr'),textoidl('H_2, hr')]
outfile = prefix + 'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536gSF'

;base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.00512.dir/h516.cosmo25cmb.1536g2MBWK.00512.halo.1', $
;        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/steps/h516.cosmo25cmb.2304g2MBWK.00512.dir/h516.cosmo25cmb.2304g2MBWK.00512.halo.1',$
;        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1']
;keys = ['1536','2304','3072']
;outfile = prefix + 'h516.cosmo25cmb.2MBK.res'



;base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.00512.dir/h516.cosmo25cmb.1536g2MBWK.00512.halo.1', $
;        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/steps/h516.cosmo25cmb.2304g2MBWK.00512.dir/h516.cosmo25cmb.2304g2MBWK.00512.halo.1',$
;        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1']
;keys = ['1536','2304','3072']
;outfile = prefix + 'h516.cosmo25cmb.2MBK.res'

base = ['h516.cosmo25cmb.768g/h516.cosmo25cmb.768g14HBWK/steps/h516.cosmo25cmb.768g14HBWK.00512.dir/h516.cosmo25cmb.768g14HBWK.00512.halo.1',$
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.halo.1',$
        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/h516.cosmo25cmb.2304g14HBWK.00512.halo.1',$
        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
;keys = ['1536','2304','3072']
keys = ['768','1536','2304','3072']
outfile = prefix + 'h516.cosmo25cmb.14HWK.res'

;base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g2MBWK/steps/h516.cosmo25cmb.1536g2MBWK.00512.dir/h516.cosmo25cmb.1536g2MBWK.00512.halo.1', $
; ;       'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g2MBWK/steps/h516.cosmo25cmb.2304g2MBWK.00512.dir/h516.cosmo25cmb.2304g2MBWK.00512.halo.1',$
;        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir/h516.cosmo25cmb.1536g14HBWK.00512.halo.1',$
; ;       'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir/h516.cosmo25cmb.2304g14HBWK.00512.halo.1',$
;        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
;keys = ['1536','2304','3072','H2 1536','H2 2304','H2 3072']
;keys = ['1536','3072','H2 1536','H2 3072']
;outfile = prefix + 'h516.cosmo25cmb.14HWK.res'

;base = ['h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g16HBWK/steps/h516.cosmo25cmb.1536g16HBWK.00512.dir/h516.cosmo25cmb.1536g16HBWK.00512.halo.1',$
;        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g16HBWK/steps/h516.cosmo25cmb.2304g16HBWK.00512.dir/h516.cosmo25cmb.2304g16HBWK.00512.halo.1']
;keys = ['1536','2304']
;outfile = prefix + 'h516.cosmo25cmb.16HWK.res'

;base = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
;        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g16HBWK/steps/h516.cosmo25cmb.2304g16HBWK.00512.dir/h516.cosmo25cmb.2304g16HBWK.00512.halo.1',$
;        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
;keys = ['no '+textoidl('H_2'), 'no '+textoidl('H_2') + ' SF',textoidl('H_2')]
;outfile = prefix + 'h516.cosmo25cmb.diffSF'

;base = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00156.dir/h516.cosmo25cmb.3072g1MBWK.00156.halo.1',$
;        'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00148.dir/h516.cosmo25cmb.3072g14HBWK.00148.halo.1',$
;        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00148.dir/h516.cosmo25cmb.2304g14HBWK.00148.halo.1',$
;        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00148.dir/h516.cosmo25cmb.1536g14HBWK.00148.halo.1']
;keys = ['no '+textoidl('H_2, hr'), textoidl('H_2, hr'), textoidl('H_2, mr'), textoidl('H_2, lr')]
;outfile = '~/h516.cosmo25cmb.diffres'

filename = prefix + base
;rotcurve, [filename], msol_per_sysmass, kpc_per_syslength, keys = keys, /column_read,/color , outfile = outfile
;stop

;-------------------------------- Many Steps
;base0 = "h516.cosmo25cmb.1536g3HBWK"
;base1 = "h516.cosmo25cmb.1536g6MbwK"
;base2 = "h516.cosmo25cmb.1536g6HBWK.jeans.prev"
;base3 = "h516.cosmo25cmb.1536g3HBWK" 
;step0 = ['00195','00240','00324','00384','00408','00456','00480','00512']
;step1 = ['00192','00240','00328','00384','00406','00455','00480','00512'];192
;step2 = ['00144','00192','00276','00336','00360','00408','00432','00464'];144
;step3 = ['00192','00240','00324','00384','00408','00456','00480','00512']
;keys = ['no H2','H2','H2 SF','MC SF']
;colors = [30,90,140,240]
prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

base0 = "h516.cosmo25cmb.1536g14HBWK"
step0 = ['00024','00048','00072','00096','00120','00144','00168','00192','00216','00240','00264','00288','00312','00336','00360','00384','00408','00432','00456','00480','00504','00512']
keys = ['new H2, lr']
colors= [240]
prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
outfile = '/astro/store/student-scratch1/christensen/MolecH/results/' + 'h516.cosmo25cmb.1536g14HBWK'

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
base0  = 'h516.cosmo25cmb.1536g14HK'
step0  = '00512'
keys   = ['no BW']
colors = [240]

cd,prefix
nsteps = N_ELEMENTS(step0) - 1
start = 0
dt = 1
msol_per_sysmass  = (fltarr(N_ELEMENTS(keys))+ 1)*2.310e15
kpc_per_syslength = (fltarr(N_ELEMENTS(keys))+ 1)*25000.

FOR i = start/dt, nsteps DO BEGIN 
;    step = i*dt
;    if (step lt 10) THEN step = '0000'+STRTRIM(step,2) ELSE BEGIN
;        if (step lt 100) THEN step = '000'+STRTRIM(step,2) ELSE step = '00'+STRTRIM(step,2)
;    ENDELSE
    filename0 = prefix + base0 + "/steps/"+base0+"."+step0[i]+".dir/" +  base0 + "." + step0[i]+'.halo.1'
;    filename1 = prefix + base1 + "/steps/"+base1+"."+step1[i]+".dir/"+ base1 + "." + step1[i]+'.halo.1'
;    filename2 = prefix + "h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/"+base2+"."+step2[i]+".dir/" + base2 + "." + step2[i]+'.halo.1' ;+STRTRIM(halo,2)
;    filename3 = prefix + base3 + "/steps_noH2SF/"+base3+"_noH2SF."+step3[i]+".dir/"+base0 + "_noH2SF." + step3[i]+'.halo.1'
    filename = [filename0];[filename1,filename3,filename0,filename2]
    stop
    rotcurve,filename,msol_per_sysmass,kpc_per_syslength,keys = keys,color=colors;,outfile = outfile+step0[i];,redshift = redshift[i];,outfile = outfile,/color
ENDFOR
end

pro rotCurve_h603
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g2HBWK/"
prefix = "/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g/"

halo = 1
msol_per_sysmass       = 1.84793e16
kpc_per_syslength      = 50000.

base = ["steps/h603.cosmo50cmb.1536g1BWK.00512.dir/h603.cosmo50cmb.1536g1BWK.00512.halo.1", $
        "h603.cosmo50cmb.2304g2bwK/h603.cosmo50cmb.2304g2bwK.00512.halo.1", $
        "h603.cosmo50cmb.2304g5bwK/h603.cosmo50cmb.2304g5bwK.00512.halo.1"]
;base = ["steps/h603.cosmo50cmb.1536g1BWK.00512.dir/h603.cosmo50cmb.1536g1BWK.00512.halo.1"]
base = ["h603.cosmo50cmb.1536g3HBWK/steps/h603.cosmo50cmb.1536g3HBWK.00120.dir/h603.cosmo50cmb.1536g3BWK.00120.halo.1"]
base = ["h603.cosmo50cmb.2304g3HBWK_00504/steps_ssource.J/h603.cosmo50cmb.2304g3HBWK.00504.00002.dir/h603.cosmo50cmb.2304g3HBWK.00504.00002.halo.1","h603.cosmo50cmb.2304g3HBWK_00504/steps_ssource/h603.cosmo50cmb.2304g3HBWK.00504.00002.dir/h603.cosmo50cmb.2304g3HBWK.00504.00002.halo.1"]
base = ['h603.cosmo50cmb.2304g3HBWK_00504/ssource.testUV/h603.cosmo50cmb.2304g3HBWK.00504.00001.dir/h603.cosmo50cmb.2304g3HBWK.00504.00001.halo.1']
base = ['h603.cosmo50cmb.2304g14HBWK/steps/h603.cosmo50cmb.2304g14HBWK.00060.dir/h603.cosmo50cmb.2304g14HBWK.00060.halo.1']

keys = ["H2","no H2, low threshold","no H2, high threshold"]
keys = ["H2"]
keys = ["Jean Crit","No Jeans Crit"]
keys = ["test UV"]
keys = ["H2"]

outfile = "/astro/net/scratch2/christensen/MolecH/results/h603g452"
outfile = "/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g3HBWK/steps/h603.cosmo50cmb.1536g3HBWK.00120.dir/h603.1536g3HBWK.00120"
outfile = "/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.2304g3HBWK_00504/h603.cosmo50cmb.2304g3HBWK.00504.00002.jeans.nojeans"
outfile = prefix+base

prefix = "/astro/store/student-scratch1/christensen/MolecH/Cosmo/"
base = ["h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.00444.dir/h603.cosmo50cmb.3072gs1MbwK.00444.halo.1",$
        "h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00444.dir/h603.cosmo50cmb.3072g14HBWK.00444.halo.1"]
keys = ["no H2",$
        "H2"]
outfile = "/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/plots/h603.cosmo50cmb.3072g.444"


filename = prefix + base
rotcurve,filename,msol_per_sysmass,kpc_per_syslength,keys = keys,/color,outfile = outfile;,/color

stop
base = "h603.cosmo50cmb.1536g3HBWK"
keys = ['h603']
steps = [37,46,48,60,61,72,84,96,108,120,128,132,168,180,192,204,216,227,228,240,252,264,271,276,288,300]
redshift = [5.895,4.971,4.805,4.005,3.950,3.431,2.995,2.651,2.370,2.136,2.000,1.937,1.482,1.363,1.256,1.159,1.071,0.997,0.990,0.915,0.846,0.782,0.746,0.721,0.665,0.612]
;steps = [300]
;redshift = [0.599]
last = 408; 168;512
dt = 1
start = 408; 168;512 
start = 0
halo = 1

nsteps = last/dt
nsteps = N_ELEMENTS(steps) - 1
cd,prefix
FOR i = start/dt, nsteps DO BEGIN 
    step = i*dt
    step = steps[i]
    if (step lt 10) THEN step = '0000'+STRTRIM(step,2) ELSE BEGIN
        if (step lt 100) THEN step = '000'+STRTRIM(step,2) ELSE step = '00'+STRTRIM(step,2)
    ENDELSE
    filename = base+"/steps/"+base+"."+step+".dir/"+base + "." + step+'.halo.'+STRTRIM(halo,2)
    print,filename
    rotcurve,[filename],msol_per_sysmass,kpc_per_syslength,keys =keys,/color,redshift = redshift[i];,outfile = outfile,/color
    stop
ENDFOR
end

pro rotCurve_twins
prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

base = ['h603.cosmo50cmb.2304g/h603.cosmo50cmb.2304g14HBWK/steps/h603.cosmo50cmb.2304g14HBWK.00240.dir/h603.cosmo50cmb.2304g14HBWK.00240.halo.1', $
        'h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00240.dir/h516.cosmo25cmb.2304g14HBWK.00240.halo.1']
keys = ["h603","h516"]
outfile = "/astro/net/scratch2/christensen/MolecH/results/h605_516.cosmo25cmb.2304g14HBWK.00240"

msol_per_sysmass = [1.84793e16,2.310e15 ]
kpc_per_syslength = [50000.,25000]

filename = prefix + base
rotcurve,[filename],msol_per_sysmass,kpc_per_syslength,keys = keys,/color;,outfile = outfile
end

pro rotCurve_res
prefix = "/astro/net/scratch2/christensen/MolecH/11M/"
;prefix = "/astro/net/nbody1/christensen/MolecH/MWHR/"

base = [  $
         "Disk_Iso_1e4/largeStar/MW_disk",$
         "Disk_Iso_1e5/largeStar/MW_disk" $
       ]
;base = [ $
;       "gas_merger0.1_split_s16n10/steps/gas_merger0.1_split_s16n10.00060.dir/gas_merger0.1_split_s16n10", $
;       "gas_merger0.1_single/steps/gas_merger0.1_single.00060.dir/gas_merger0.1_single" $
;]

keys = ['1e4 GP','1e5 GP' ];k
;keys = ['Split','Single' ]

outfile = "/astro/net/scratch2/christensen/MolecH/results/11Mdisk_res_"
; outfile = "/astro/net/scratch2/christensen/MolecH/results/gas_merger_res"

msol_per_sysmass         = 1.36e17
kpc_per_syslength        = 1e5
;msol_per_sysmass         = 2.3262e5
;kpc_per_syslength        = 1.0

last = 3;512
dt = 1
start = 3;512 
;last = 60;512
;dt = 1
;start = 60;512

FOR i = start/dt, last/dt DO BEGIN 
    if (i*dt lt 10) THEN step = '0000'+STRTRIM(i*dt,2) ELSE BEGIN
        if (i*dt lt 100) THEN step = '000'+STRTRIM(i*dt,2) ELSE step = '00'+STRTRIM(dt*i,2)
    ENDELSE
    filename = prefix+base+"."+step
    rotcurve,[filename],msol_per_sysmass,kpc_per_syslength,keys = keys,outfile = outfile + step
    stop
ENDFOR
END



pro rotCurve, files, dMsolUnit, dKpcUnit, column_read=column_read, cosmo = cosmo, outfile = outfile, keys = keys, color = color, redshift = redshift
;*********************************** Units ****************************************
hubble = 73.0
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
dDelta = 0.003
sec_per_year = 31556926
molec_weight = (0.76*1 + 0.24*4.0)
show_den = 1
show_temp = 1
show_phase = 0
show_frac = 1
show_rot = 1
show_line = 1
show_prof = 1
show_stellarprof = 1
show_sfr = 1
show_column = 0
show_column_check = 0
SHOW_LWFIELD = 0
dtime = 100.0e6
if not KEYWORD_SET(column_read) THEN column_read = 0 

;************************************ Plot Format ************************************
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outfile) THEN BEGIN
    set_plot,'ps' 
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
;    !p.font=0 
    cb_charsize = 0.75
    l_charsize = 0.75
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    fgcolor = 0
ENDIF ELSE BEGIN
    set_plot,'x'
    nbins=100.0
    linestyles = [0,2]
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
    cb_charsize = 1.0
    l_charsize = 1.0
    !X.MARGIN = [12,3]
    !Y.MARGIN = [6,2]
    fgcolor = 255
ENDELSE
!p.multi = 0

IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    if color[0] eq 1 then  colors = (findgen(N_ELEMENTS(files)) + 1)*240/N_ELEMENTS(files) else colors = color
    thicks = fltarr(N_ELEMENTS(files)) + 2
    symbols = fltarr(N_ELEMENTS(files)) + 4
ENDIF ELSE BEGIN
    loadct,0    
    colors =  fltarr(N_ELEMENTS(files)) + fgcolor;(findgen(N_ELEMENTS(files)) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(files)) + 5
    thicks = (findgen(N_ELEMENTS(files)) + 1)*6/N_ELEMENTS(files) - 1
    symbols = (findgen(N_ELEMENTS(files))+1)*2
ENDELSE

IF NOT (KEYWORD_SET(legend_t)) THEN legend_t = ''
names = files
halo_stats = REPLICATE({halostats, tmass:1.2e10, dmass:1.2e10, smass:1.2e10, gmass:1.2e10, HImass:1.2e10, H2mass:1.2e10, gasSD:1.2e10, starSD:1.2e10, rad:1.2e10, v_f:1.2e10}, N_ELEMENTS(files))

;*********************************** Loading Sims ***********************************
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    filename = files[ifile]
    names[ifile] = (STRSPLIT(filename,'/',/extract))[N_ELEMENTS(STRSPLIT(filename,'/')) - 1]
;   IF NOT (KEYWORD_SET(redshift)) THEN BEGIN;redshift  = 0
;       ind = STRSPLIT(filename,'/')
;       base = STRMID(filename,0,ind[N_ELEMENTS(ind) - 1])
;       spawn,'ls '+base+'*.AHF_halos',result
;       segments = strsplit(result,'.',/extract)
;       redshift = segments[N_ELEMENTS(segments) - 3] + '.'  +segments[N_ELEMENTS(segments) - 2]
;       redshift = STRMID(redshift,1,STRLEN(redshift) - 1)
;   ENDIF
 ;  IF (N_ELEMENTS(redshift) GT 1) THEN a = 1.0/(1.0 + redshift[ifile]) ELSE a = 1.0/(1.0 + redshift[0])
    IF (N_ELEMENTS(dKpcUnit) GT 1) THEN kpc_per_syslength = dKpcUnit[ifile] ELSE kpc_per_syslength = (dKpcUnit)[0]
    IF (N_ELEMENTS(dMsolUnit) GT 1) THEN msol_per_sysmass = dMsolUnit[ifile] ELSE msol_per_sysmass = dMsolUnit[0]
    timeunit = SQRT((cm_per_kpc*kpc_per_syslength)^3/(gm_per_msol*msol_per_sysmass)/6.67d-8)/sec_per_year
    IF KEYWORD_SET(cosmo) THEN BEGIN
        filenametip = filename + '.halo.'+STRTRIM(cosmo,2)+'.std'
        grpfile = filename+'.amiga.grp'
        grp = read_lon_array(grpfile)
        IF (FILE_TEST(filenametip)) THEN BEGIN  ;If the halo has already been selected (saves lots of time)
            rtipsy,filename,h,g,d,s,/justhead
            a = h.time
            ind = where(grp eq cosmo,comp=indcomp)
            indg = ind[where(ind lt h.ngas)]
            rtipsy,filenametip,hhalo,g,d,s 
        ENDIF ELSE BEGIN  ;If the halo hasn't been selected, then select and rotate it
            rtipsy,filename,h,g,d,s
            a = h.time
            ind = where(grp eq cosmo,comp=indcomp)
            inds = ind[where(ind ge h.ngas+h.ndark)]-h.ngas-h.ndark
            indg = ind[where(ind lt h.ngas)]
            indd = ind[where(ind ge h.ngas and ind lt h.ngas+h.ndark)]-h.ngas

            s = s[inds]
            g = g[indg]
            d = d[indd]
            rtipsy, filename+'.amiga.gtp', h1,g1,d1,s1
        
            cx= s1[0].x
            cy= s1[0].y
            cz= s1[0].z
        
            g.x = g.x - cx
            g.y = g.y - cy
            g.z = g.z - cz
            d.x = d.x - cx
            d.y = d.y - cy
            d.z = d.z - cz
            s.x =  s.x - cx
            s.y =  s.y - cy
            s.z =  s.z - cz
            dist_units=kpc_per_syslength*h.time ; multiplies the distance units, which you fed it, by the expansion         
            limit=5*h.time/dist_units ;use stars within this limit (physical kpc/dist_units)                    
            align,g,d,s,limit
        ENDELSE
        readarr,filename+".HI",h,HI_prime,/ascii
        IF (FILE_TEST(filename+".H2")) THEN  readarr,filename+".H2",h,H2_prime,/ascii else H2_prime = fltarr(N_ELEMENTS(HI_prime))
        IF(show_column_check OR show_column) THEN BEGIN
            IF (FILE_TEST(filename+'.smoothlength')) THEN BEGIN
                readarr,filename+'.smoothlength',h,smoothlengths_L,/ascii
                smoothlengths = (smoothlengths_L*kpc_per_syslength*cm_per_kpc)
            ENDIF ELSE BEGIN
                show_column = 1
                show_column_check = 0
            ENDELSE
            IF KEYWORD_SET(column_read && FILE_TEST(filename+'.correL')) then BEGIN
;                readarr,filename+'.column',h,column_L,/ascii
                readarr,filename+'.correL',h,column_L,/ascii
                column = column_L[0:h.ngas - 1] 
                column = column*cm_per_kpc*kpc_per_syslength
                mach = smoothlengths/column
            ENDIF ELSE BEGIN
                IF (FILE_TEST(filename+'.shear')) THEN BEGIN
                    readarr,filename+'.shear',h,mach,/ascii
                    column = (1.0/mach*kpc_per_syslength*cm_per_kpc)
                ENDIF ELSE BEGIN
;                    show_column = 0
                    show_column_check = 0
                ENDELSE
            ENDELSE
    length = smoothlengths ;column 
        ENDIF
    kpc_per_syslength = kpc_per_syslength*a        
    ENDIF ELSE BEGIN
        rtipsy,filename,h,g,d,s
        readarr,filename+".HI",h,HI_prime,/ascii
        IF (FILE_TEST(filename+".H2")) THEN  readarr,filename+".H2",h,H2_prime,/ascii else H2_prime = fltarr(N_ELEMENTS(HI_prime))
        IF (FILE_TEST(filename+'.lw')) THEN BEGIN
            readarr,filename+".lw",h,lw_prime,/ascii
            lw = lw_prime[0:h.ngas - 1]*1d30/cm_per_kpc^2/kpc_per_syslength^2/h.time/h.time
        ENDIF
        HI_frac = HI_prime[0:h.ngas - 1]
        H2_frac = H2_prime[0:h.ngas - 1]
        coolon = g.dens*0 + 1.0
        IF (FILE_TEST(filename+".coolontime")) THEN  BEGIN
            readarr,filename+".coolontime",h,coolontime,/ascii,part = 'gas'
            coolontime = coolontime*TIMEUNIT
            snaptime = MAX(s.tform)*TIMEUNIT
            coolon[where(coolontime - snaptime gt 0)] = 0
        ENDIF ELSE coolon = g.dens*0 + 1.0
        IF(show_column_check OR show_column) THEN BEGIN
            IF (FILE_TEST(filename+'.smoothlength')) THEN BEGIN
                readarr,filename+".smoothlength",h,smoothlengths_L,/ascii
                smoothlengths = smoothlengths_L[0:h.ngas - 1]*kpc_per_syslength*cm_per_kpc
                column = smoothlengths
            ENDIF ELSE BEGIN
                show_column = 1
                show_column_check = 0
            ENDELSE 
            IF KEYWORD_SET(column_read && FILE_TEST(filename+'.correL')) then BEGIN
;                readarr,filename+".column",h,column_L,/ascii
                readarr,filename+".correL",h,column_L,/ascii
                column = column_L[0:h.ngas - 1] 
                column = column*cm_per_kpc*kpc_per_syslength
                mach = smoothlengths/column
            ENDIF ELSE BEGIN
                IF (FILE_TEST(filename+'.shear')) THEN BEGIN
                    readarr,filename+".shear",h,mach,/ascii
                    column = 1.0/mach[0:h.ngas - 1]*kpc_per_syslength*cm_per_kpc
                ENDIF ELSE BEGIN
                    show_column = 1
                    show_column_check = 0
                ENDELSE
            ENDELSE    
            length = smoothlengths ;column
        ENDIF
        print,"Data read"
    ENDELSE
    dens_convert =  msol_per_sysmass * gm_per_msol * 5.9753790e+23/kpc_per_syslength^3/cm_per_kpc^3
    vel_convert = hubble * (kpc_per_syslength / 1000.0)/2.894405
    total_H = fltarr(N_ELEMENTS(HI_frac))
    total_He = fltarr(N_ELEMENTS(HI_frac))
    total_He[where(g.zmetal le 0.1, COMPLEMENT = comple)] = (0.236 + 2.1*g[where(g.zmetal le 0.1)].zmetal)/4.0
    IF(comple[0] ne -1) THEN total_He[comple] = (-0.446*(g[comple].zmetal - 0.1)/0.9 + 0.446)/4.0 ;
    total_H = 1.0 - total_He*4.0 - g.zmetal 
    HI = HI_frac*g.mass*msol_per_sysmass
    H2 = H2_frac*g.mass*msol_per_sysmass
    rho = alog10(g.dens * dens_convert * total_H) ;H Atoms per CC 
    t = alog10(g.tempg)
    v = g.vx*vel_convert
    r = sqrt(g.x*g.x + g.y*g.y)*kpc_per_syslength
    halo_stats[iFile].tmass = (TOTAL(g.mass) + TOTAL(d.mass) + TOTAL(s.mass))*msol_per_sysmass
    halo_stats[iFile].dmass = TOTAL(d.mass)*msol_per_sysmass
    halo_stats[iFile].gmass = TOTAL(g.mass)*msol_per_sysmass
    halo_stats[iFile].smass = TOTAL(s.mass)*msol_per_sysmass
    halo_stats[iFile].HImass = TOTAL(HI)
    halo_stats[iFile].H2mass = TOTAL(H2)
;********************************* Column Check***********************************************************
    !p.multi = 0
    IF (show_column_check) THEN BEGIN
        IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+STRTRIM(ifile,2)+'_smooth.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  ELSE window,0,xsize = 712,ysize = 392
        dens_ind = where(rho gt 0)
        smoothl_hist = histogram(alog10(smoothlengths),locations = x,min = 14,max = 24,nbins = 50)
        column_hit = histogram(alog10(column),locations = x,min = 14,max = 24,nbins = 50)
        columnden_hit = histogram(alog10(g.dens*dens_convert*total_H*column),locations = x,min = 14,max = 24,nbins = 50)
        IF (dens_ind[0] ne -1) THEN BEGIN
            smoothl_hit_dens = histogram(alog10(smoothlengths[dens_ind]),locations = x,min = 14,max = 24,nbins = 50)
            column_hit_dens = histogram(alog10(column[dens_ind]),locations = x,min = 14,max = 24,nbins = 50)
        ENDIF
        plot,x,smoothl_hist,psym = 10,xtitle = "Log(h) [cm]",title = names[ifile],xrange = [16,24],thick = 0.5
        oplot,x,column_hit,linestyle = 2,psym = 10,thick = 0.5
        oplot,x,columnden_hit,linestyle = 1,psym = 10
        IF (dens_ind[0] ne -1) THEN BEGIN
            oplot,x,smoothl_hit_dens,linestyle = 0,psym = 10,thick = 3
            oplot,x,column_hit_dens,linestyle = 2,psym = 10,thick = 3
            legend,["Smoothing Length","Column Length","Column Density","High Density Smooth Lengths","High Density Column Lengths"],linestyle=[0,2,1,0,2],charsize = l_charsize,thick = [1,1,1,3,3]
        ENDIF ELSE legend,["Smoothing Length","Column Length","Column Density"],linestyle=[0,2,1],charsize = l_charsize
        IF (KEYWORD_SET(outfile)) THEN device,/close 
    ENDIF

;********************************* H2 Fraction Column ***********************************************************
    IF (show_column AND MAX(H2_frac) ne 0) THEN BEGIN
        readcol,'/astro/users/christensen/code/MolecH/Wolfire08.dat',starw,namew,nhw,nh_erw,nh2w,nh2_erw,ncIw,ncI_erw,ncIIw,ncII_erw,avw,logfH2w,logfcIw,refw,format='A10,A8,F,F,F,F,F,F,F,F,F,F,F,I'
        readcol,'/astro/users/christensen/code/MolecH/Gillmon06.dat',nameg,nh2g,nhg,refg,logfH2g,T01,Texc,format = 'A11,F,F,I,F,I,I,I'
        IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+STRTRIM(ifile,2)+'_frac_columnMass.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  ELSE window,1,xsize = 712,ysize = 392
        plot,alog10(g.dens*dens_convert*total_H*length),H2_frac/(total_H/2.0), xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),title = names[ifile],psym=3,xstyle=1,ystyle=1,symsize =0.5,xrange = [19,24],yrange = [1e-6,1.0],/ylog ;,title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
        
;        contour_plus,alog10(g.dens*dens_convert*total_H*smoothlengths),H2_frac/(total_H/2.0),threshold=3000,nlevels=10,XBINSIZE=0.1,YBINSIZE=0.05,YMIN=0,YMAX=1.0, xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [0, 1.1],title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
;       oplot,[-4,4],[1,1]
;        contour_plus,alog10(g.dens*dens_convert*total_H*column),H2_frac/(total_H/2.0),threshold=300000,nlevels=10,XBINSIZE=0.1,YBINSIZE=0.05,YMIN=1e-6,YMAX=1.1,xmax = 24, xmin = 19, xtitle = textoidl("N_{HI} + 2N_{H_2} (cm^{-2})"),ytitle = textoidl('f_{H_2}'),yrange = [1e-6,1.1],title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years',/ylog

;        oplot,nhw,10^logfH2w,psym = 1,color = 240
;        oplot,nhg,10^logfH2g,psym = 4,color = 240
;        legend,['Simulated Data','FUSE, Gillmon et al. 06','Wolfire et al. 08'],psym = [3,4,1],color=[0,240,240],/right,/bottom,charsize = l_charsize


        column_den= 10^(findgen(500)/500.0*(23-19) + 19.0)
  ;     stop
;        z = 0.00158489 ;MEAN(g[where(frac_ind gt 1e-6)].zmetal)
        ZSOLAR = 0.0177;0.0130215
        z = zsolar*0.01
        sigma_d = 2d-21
        shield = exp(-1.0*sigma_d*z/ZSOLAR*(column_den))
;        oplot,alog10(column_Den),1.0-shield
;        shield = exp(-1.0*sigma_d*(column_den))
;        oplot,alog10(column_Den),1.0-shield,linestyle = 2
        lte_H2_Dat = column_den
        cl = 1e20 ;5.29e19
        if 1 then begin
        lw_lte = 1e4
        FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
          lte_H2_Dat[i] = (lte_H2(column_den[i]/cl/MEAN(total_H),z,cl,lw_lte,MEAN(total_H)))[0]
        oplot,alog10(column_Den),lte_H2_Dat/MEAN(total_H)*2.0,linestyle = 2,color =240*(alog10(lw_lte)-3)/(10-3)
        lw_lte = 1e5
        FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
          lte_H2_Dat[i] = (lte_H2(column_den[i]/cl/MEAN(total_H),z,cl,lw_lte,MEAN(total_H)))[0]
        oplot,alog10(column_Den),lte_H2_Dat/MEAN(total_H)*2.0,linestyle = 2,color =240*(alog10(lw_lte)-3)/(10-3)
        lw_lte = 1e6
        FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
          lte_H2_Dat[i] = (lte_H2(column_den[i]/cl/MEAN(total_H),z,cl,lw_lte,MEAN(total_H)))[0]
        oplot,alog10(column_Den),lte_H2_Dat/MEAN(total_H)*2.0,linestyle = 2,color =240*(alog10(lw_lte)-3)/(10-3)
        lw_lte = 1e7
        FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
          lte_H2_Dat[i] = (lte_H2(column_den[i]/cl/MEAN(total_H),z,cl,lw_lte,MEAN(total_H)))[0]
        oplot,alog10(column_Den),lte_H2_Dat/MEAN(total_H)*2.0,linestyle = 2,color =240*(alog10(lw_lte)-3)/(10-3)
        lw_lte = 1e8
        FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
          lte_H2_Dat[i] = (lte_H2(column_den[i]/cl/MEAN(total_H),z,cl,lw_lte,MEAN(total_H)))[0]
        oplot,alog10(column_Den),lte_H2_Dat/MEAN(total_H)*2.0,linestyle = 2,color =240*(alog10(lw_lte)-3)/(10-3)
        lw_lte = 1e9
        FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
          lte_H2_Dat[i] = (lte_H2(column_den[i]/cl/MEAN(total_H),z,cl,lw_lte,MEAN(total_H)))[0]
        oplot,alog10(column_Den),lte_H2_Dat/MEAN(total_H)*2.0,linestyle = 2,color =240*(alog10(lw_lte)-3)/(10-3) 
        lw_lte = 2e10
        FOR i = 0, N_ELEMENTS(column_den) - 1 DO $
          lte_H2_Dat[i] = (lte_H2(column_den[i]/cl/MEAN(total_H),z,cl,lw_lte,MEAN(total_H)))[0]
        oplot,alog10(column_Den),lte_H2_Dat/MEAN(total_H)*2.0,linestyle = 2,color =240*(alog10(lw_lte)-3)/(10-3)
        endif
        IF 1 THEN BEGIN
            ind = where(H2_frac/total_H*2 ge 5e-7 AND coolon)
            dens_ind = alog10(g[ind].dens*dens_convert*total_H[ind]*length[ind]) ;smoothlengths[ind])
            frac_ind = H2_frac[ind]/(total_H[ind]/2.0) 
            ;where(alog10(g.zmetal) lt 0 AND alog10(g.zmetal) gt -1)            

            colorarr = alog10(lw[ind])
                                ;  colorarr = alog10(g.zmetal)
 ;           colorarr = alog10(column[ind])
            min_colorarr = 4 ;MIN(colorarr)
            max_colorarr = 10 ;MAX(colorarr)
;            min_colorarr = 3;25.0;-4.5
;            max_colorarr = 9;30.0;2
            ;colorarr[where(colorarr lt min_colorarr)] = min_colorarr
            ;colorarr[where(colorarr gt max_colorarr)] = max_colorarr
            colors_pd = FIX((colorarr - min_colorarr)/(max_colorarr - min_colorarr)*254)
            if (where(colorarr lt min_colorarr))[0] ne -1 THEN colors_pd[where(colorarr lt min_colorarr)] = 255
            if (where(colorarr gt max_colorarr))[0] ne -1 THEN colors_pd[where(colorarr gt max_colorarr)] = 255
            for ii=0LL,LONG64(N_ELEMENTS(g[ind].dens)) - 1 do $
              oplot,[dens_ind[ii],dens_ind[ii]], [frac_ind[ii],frac_ind[ii]],psym = 3,color=colors_pd[ii]
            if KEYWORD_SET(outfile) THEN position = [0.88, 0.14, 0.95, 0.92] ELSE position = [0.88, 0.1, 0.95, 0.9]     
            colorbar,maxrange = max_colorarr,minrange = min_colorarr,/vertical,position = position,format = '(F10.1)'
            colorarr_hist = histogram(colorarr,locations = x,nbins = 100,min = min_colorarr,max = max_colorarr)
            plot,colorarr_hist,x,psym = 10,position = position,/NOERASE,XRANGE = [0,MAX(colorarr_hist)],YRANGE = [min_colorarr, max_colorarr],xstyle = 9,ystyle = 9,xticks = 1,yticks = 1,XTICKFORMAT='(A1)',YTICKFORMAT='(A1)'
;             IF (KEYWORD_SET(outfile)) THEN BEGIN
;                device,/close
;                device,filename=outfile+STRTRIM(ifile,2)+'frac_columnMass.eps',/color,bits_per_pixel= 8,/times 
;            ENDIF ELSE window,10
;            colorarr_hist = histogram(colorarr,locations = x,nbins = 100,min = min_colorarr,max = max_colorarr)
;            plot,x,colorarr_hist,psym = 10,xrange = [min_colorarr,max_colorarr],yrange = [0,N_ELEMENTS(ind)*0.05]
        ENDIF
        IF (KEYWORD_SET(outfile)) THEN BEGIN
            device,/close
        ENDIF
    ENDIF

    IF(show_LWfield) THEN BEGIN
        window,3
        plot,g.dens*dens_convert*total_H,lw/1e6,psym = 3,/ylog,/xlog,xrange = [0.1,1e3],yrange = [0.001,10000],ytitle = 'f_1000A/f_Draine',xtitle = 'n_H (cm^-3)' 
    ENDIF

;********************************* Phase Diagram ***********************************************************
    IF (show_phase AND (show_frac AND ((WHERE(H2 ne 0))[0] ne -1))) THEN BEGIN
        IF (KEYWORD_SET(outfile)) THEN BEGIN
            device,filename=outfile+STRTRIM(ifile,2)+'_phase_frac.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  
            position = [0.16, 0.90, 0.45, 0.93]
        ENDIF ELSE BEGIN
            window,2,xsize = 712,ysize = 392
            position = [0.10, 0.88, 0.42, 0.95]
        ENDELSE
        !p.multi = [0,2,1]
    ENDIF ELSE BEGIN
        IF (show_phase) THEN BEGIN
            IF (KEYWORD_SET(outfile)) THEN BEGIN
                device,filename=outfile+STRTRIM(ifile,2)+'_phase.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
                position = [0.10, 0.88, 0.90, 0.95]
            ENDIF ELSE BEGIN
                window,2 ,xsize = 712,ysize = 392
                position = [0.10, 0.88, 0.90, 0.95]
            ENDELSE
        ENDIF
        IF (show_frac AND ((WHERE(H2 ne 0))[0] ne -1)) THEN BEGIN
            IF (KEYWORD_SET(outfile)) THEN  device,filename=outfile+STRTRIM(ifile,2)+'_frac_.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2  ELSE window,2,xsize = 712,ysize = 392
        ENDIF
    ENDELSE
    
    
    IF (show_phase) THEN BEGIN
;        IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+STRTRIM(ifile,2)+'phase.eps',/color,bits_per_pixel= 8,/times ELSE window,1
        logfrac = alog10(H2_frac/(total_H/2.0))
        minlogfrac = -6         ;MIN(alog10(H2_frac/(total_H/2.0)))
        maxlogfrac = 0 ;-0.0518840 ;MAX(alog10(H2_frac/(total_H/2.0))) 
        colors_pd = FIX((logfrac - minlogfrac)/(maxlogfrac - minlogfrac)*254)
        IF KEYWORD_SET(outfile) THEN colors_pd[where(logfrac lt minlogfrac)] = 1 ELSE colors_pd[where(logfrac lt minlogfrac)] = 255
        plot,g.dens*dens_convert*total_H,g.tempg,psym = 3,xtitle = textoidl("n_H (cm^{-3})"),ytitle = "T (K)",xrange = [1e-6,10000],yrange = [1e1, 1e7],/xlog,/ylog,xstyle = 1;,title = names[ifile] ;,title='Phase Diagram '+legend_t+' -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years'
        for ii=0LL,LONG64(N_ELEMENTS(t)) - 1 do oplot,[g[ii].dens*dens_convert*total_H[ii],g[ii].dens*dens_convert*total_H[ii]],[g[ii].tempg,g[ii].tempg],psym = 3,color=colors_pd[ii]
;   contour_plus,alog10(g.dens*dens_convert),alog10(g.tempg),weight = H2_frac/(total_H/2.0),threshold = 0.1,xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = "LOG(T) (K)",xrange = [-4,4.5],yrange = [1, 7],title='Phase Diagram -- '+strtrim(step),levels = levels
        
        colorbar,maxrange = maxlogfrac,minrange = minlogfrac,format = '(F10.1)',position = position,charsize = cb_charsize
;        IF (KEYWORD_SET(outfile)) THEN BEGIN
;            device,/close
;        ENDIF
    ENDIF
    
;********************************* H2 Fraction ***********************************************************
    IF (show_frac AND ((WHERE(H2 ne 0))[0] ne -1)) THEN BEGIN
;       IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+STRTRIM(ifile,2)+'frac_.eps',/color,bits_per_pixel= 8,/times ELSE window,2
        plot,alog10(g.dens*dens_convert*total_H),H2_frac/(total_H/2.0),xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = textoidl('f_{H_2}'),xrange = [0,4],yrange = [0, 1.1],psym=4,xstyle=1,ystyle=1 ;,title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
;       contour_plus,alog10(g.dens*dens_convert*total_H),H2_frac/(total_H/2.0),threshold=300000,nlevels=10,XBINSIZE=0.1,YBINSIZE=0.05,YMIN=0,YMAX=1.0, xtitle = textoidl("LOG(n_H) (cm^{-3})"),ytitle = textoidl('f_{H_2}'),xrange = [-1,4],yrange = [0, 1.1],title='Fraction H2 -- '+ strtrim(i*dt*dDelta*timeunit,2) + ' years'
        
;       ind_high = where(column le 1e19)
;       ind_mid = where(column gt 1e19 AND column le 1e20)
;       ind_low = where(column gt 1e20)        
;       oplot,alog10(g[ind_high].dens*dens_convert*total_H[ind_high]),H2_frac[ind_high]/total_H[ind_high]*2.,color = 240,psym = 3
;       oplot,alog10(g[ind_mid].dens*dens_convert*total_H[ind_mid]),H2_frac[ind_mid]/total_H[ind_mid]*2.,color = 100,psym = 3
;       oplot,alog10(g[ind_low].dens*dens_convert*total_H[ind_low]),H2_frac[ind_low]/total_H[ind_low]*2.,color = 60,psym = 3
        oplot,[-4,4],[1,1]
;       stop
;       IF (KEYWORD_SET(outfile)) THEN device,/close
    ENDIF
    IF (KEYWORD_SET(outfile) AND (show_phase OR (show_frac AND ((WHERE(H2 ne 0))[0] ne -1)))) THEN device,/close
    !p.multi = 0
    maxr = 12. ;8.

;********************************* Rotation Curve **********************************************************
    IF (show_rot) THEN BEGIN
        maxdistance = 8
        bind = maxdistance/kpc_per_syslength/nbins
        vc_x = (findgen(nbins)+1)*bind
        tmass = fltarr(nbins)
        density = fltarr(nbins)

        distancesg = SQRT(g.x^2.0 + g.y^2.0 + g.z^2.0)
        IF N_ELEMENTS(s) gt 0 THEN distancess = SQRT(s.x^2.0 + s.y^2.0 + s.z^2.0)
        distancesd = SQRT(d.x^2.0 + d.y^2.0 + d.z^2.0)
        IF N_ELEMENTS(s) gt 0 THEN distances = [distancesg,distancess,distancesd] ELSE distances = [distancesg,distancesd]
        
        massg = g.mass*msol_per_sysmass
        IF N_ELEMENTS(s) gt 0 THEN masss = s.mass*msol_per_sysmass
        massd = d.mass*msol_per_sysmass
        IF N_ELEMENTS(s) gt 0 THEN mass = [massg,masss,massd] ELSE mass = [massg,massd]

        grav = 6.67e-8
        FOR i = 0, nbins-1 DO BEGIN
            ind = where(distances LE vc_x[i])
            if (ind[0] ne -1)then tmass[i] = TOTAL(mass[ind])else tmass[i] = 0
        ENDFOR
        velocitycurve = SQRT(tmass*gm_per_msol*grav/vc_x/kpc_per_syslength/cm_per_kpc)/1e5
        vc_x = vc_x*kpc_per_syslength
        IF iFile eq 0 THEN vc_y = fltarr(N_ELEMENTS(files),N_ELEMENTS(vc_x))
        vc_y[iFile,*] = velocitycurve
        halo_stats[iFile].v_f = MAX(vc_y[*,where(vc_x gt 3)])
    ENDIF

;********************************* Line Profile ***********************************************************
    IF (show_line) THEN BEGIN
        minv = -300.
        maxv = 300.
        range = maxv - minv
        v_bin = findgen(nbins)*range/(nbins) - range/2.0
        y_H2 = weighted_HISTOGRAM(v,weight=H2,binsize=range/nbins,min=minv,max=maxv)
        y_HI = weighted_HISTOGRAM(v,weight=HI,binsize=range/nbins,min=minv,max=maxv)
        IF iFile eq 0 THEN v_y = fltarr(2,N_ELEMENTS(files),N_ELEMENTS(y_H2))
        v_y[0,iFile,*] = y_H2
        v_y[1,iFile,*] = y_HI   
    ENDIF

;********************************* Radial Profile ***********************************************************
    IF (show_prof) THEN BEGIN
        minr = 0
        range = maxr - minr
        r_bin = findgen(nbins)*range/(nbins)
        y_H2 = weighted_HISTOGRAM(r,weight=H2,binsize=range/nbins,min=minr,max=maxr)
        y_HI = weighted_HISTOGRAM(r,weight=HI,binsize=range/nbins,min=minr,max=maxr)
        dr = range/(nbins)
        area = ((r_bin + dr)*(r_bin + dr) - r_bin*r_bin)*1e6
        IF iFile eq 0 THEN prof_y = fltarr(2,N_ELEMENTS(files),N_ELEMENTS(y_H2))
        prof_y[0,iFile,*] = y_H2/area
        prof_y[1,iFile,*] = y_HI/area       
;        stop
    ENDIF

;********************************* Temperature Distro ***********************************************************
    IF (show_temp) THEN BEGIN
        mint= 1
        maxt = 4.5
        range = maxt - mint
        t_bin = findgen(nbins)*range/(nbins) + mint
        y_H2 = weighted_HISTOGRAM(t,weight=H2,binsize=range/nbins,min=mint,max=maxt)
        y_HI = weighted_HISTOGRAM(t,weight=HI,binsize=range/nbins,min=mint,max=maxt)
        IF iFile eq 0 THEN t_y = fltarr(2,N_ELEMENTS(files),N_ELEMENTS(y_H2))
        t_y[0,iFile,*] = y_H2
        t_y[1,iFile,*] = y_HI   
    ENDIF
 
;********************************* Density Distro ***********************************************************
    IF (show_den) THEN BEGIN    
        minrho = -5
        maxrho = 4
        range = maxrho - minrho
        rho_bin = findgen(nbins)*range/(nbins) + minrho
        y_H2 = weighted_HISTOGRAM(rho,weight=H2,binsize=range/nbins,min=minrho,max=maxrho)
        y_HI = weighted_HISTOGRAM(rho,weight=HI,binsize=range/nbins,min=minrho,max=maxrho)
        IF iFile eq 0 THEN rho_y = fltarr(2,N_ELEMENTS(files),N_ELEMENTS(y_H2))
        rho_y[0,iFile,*] = y_H2
        rho_y[1,iFile,*] = y_HI   
    ENDIF

;************************************ Stellar Prof **********************************
    IF (show_stellarprof) THEN BEGIN
        minr = 0
        range = maxr - minr
        star  = s
        star.mass = s.mass*msol_per_sysmass
        star.x = s.x*kpc_per_syslength
        star.y = s.y*kpc_per_syslength
        star.z = s.z*kpc_per_syslength        
        prof = prof(star,'star',MAX(s.tform),nbins = nbins, rmax = maxr, rmin = minr)
        IF iFile eq 0 THEN prof_stellar = fltarr(2,N_ELEMENTS(files),nbins)
        prof_stellar[0,iFile,*] = prof.rbins
        prof_stellar[1,iFile,*] = prof.rho 

        tcurrent = max(s.tform*timeunit)
        massform = MAX(s.mass)*msol_per_sysmass
        fitpar =  [0,0,0,0,0]
        ind = where(prof.rbins gt 1 and prof.rbins lt 8)   
        dblexpfit, prof.rbins[ind], prof.rho[ind], $
          prof.rho[ind]/sqrt(prof.num_in_bin[ind]), fitpar, red_chisq=exp2chi, /no_guess
        tclip = tcurrent - dtime
        halo_stats[iFile].rad = 3;fitpar[2]
        rform=SQRT(s.x^2+s.y^2)*kpc_per_syslength
        rgas=SQRT(g.x^2+g.y^2)*kpc_per_syslength
        indform=WHERE(rform LT halo_stats[iFile].rad AND s.tform*timeunit GT tclip,nindform); AND ABS(s.z)*kpcunit LT diskheight 
        indgas=WHERE(rgas LT halo_stats[iFile].rad,nindgas) ;AND ABS(g.z)*kpcunit LT diskheight 
        halo_stats[iFile].gasSD = TOTAL((HI[indgas]))/(!PI*halo_stats[iFile].rad^2*1.e6)
        halo_stats[iFile].starSD = TOTAL(massform*nindform)/(!PI*halo_stats[iFile].rad^2*dtime)
    ENDIF        

;********************************* Star Formation Rate ***********************************************************
    IF (show_sfr) THEN BEGIN    
        IF iFile eq 0 THEN BEGIN
            mintime = 0
            maxtime = MAX([MAX(s.tform)*timeunit/1e9,14])
            maxtime = MAX(s.tform)*timeunit/1e9
            rangesfr = maxtime - mintime
            tbinsize = 0.05
            ntbins = range/tbinsize
            sfr_bin = findgen(ntbins)*rangesfr/(ntbins)
        ENDIF
        y_s = weighted_HISTOGRAM(s.tform*timeunit/1e9,weight=massform/(rangesfr/ntbins*1e9),binsize=rangesfr/ntbins,min=mintime,max=maxtime)
        IF iFile eq 0 THEN sfr_y = fltarr(1,N_ELEMENTS(files),N_ELEMENTS(y_s))
        sfr_y[0,iFile,*] = y_s
      ENDIF

    print,'Temperature (Mean, Min, Max): ',MEAN(g.tempg),MINMAX(g.tempg)
    print,'Density (Mean, Min, Max): ',MEAN(rho),MINMAX(rho)
    print,'H2 (Mean, Min, Max): ',MEAN(H2_frac/(total_H/2.0)),MINMAX(H2_frac/(total_H/2.0))
    print,'Total Mass in H2 + HI + He: ',TOTAL(1.36*(HI_frac + 2.0*H2_frac)/(total_H)*g.mass*msol_per_sysmass)
    print,'Total Mass in H2 + HI: ',TOTAL((HI_frac + 2.0*H2_frac)/(total_H)*g.mass*msol_per_sysmass)   
    print,'Total Mass in H2: ',TOTAL(2.0*H2_frac/(total_H)*g.mass*msol_per_sysmass)
    print,'Percent Mass in Cold ISM: ',TOTAL(H2_frac/(total_H/2.0)*g.mass*msol_per_sysmass)/TOTAL(g.mass*msol_per_sysmass)*100.0
    print,'Total Stellar Mass: ',TOTAL(s.mass)*msol_per_sysmass
    print,' '
    print,'H2/HI: ',TOTAL(2.0*H2_frac/(total_H)*g.mass*msol_per_sysmass)/TOTAL(HI_frac/(total_H)*g.mass*msol_per_sysmass)
    IF NOT (KEYWORD_SET(outfile)) THEN BEGIN
        IF N_ELEMENTS(keys) EQ 1 THEN wait,0.5 ELSE stop
    ENDIF
ENDFOR


IF NOT (KEYWORD_SET(keys)) THEN keys = names
!p.multi = [0,2,3]
IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+'_rotCurve.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 25,xoffset =  2,yoffset = 2 ELSE window,3,xsize = 712,ysize = 800;392
;********************************* Rotation Curve ***********************************************************
;IF (show_rot) THEN BEGIN
    maxvc = 70; 250 ;70 ;120
IF (KEYWORD_SET(outfile)) THEN BEGIN
;    device,filename=outfile+'rotCurve_.eps',/color,bits_per_pixel=
;    8,/times 
    IF KEYWORD_SET(COLOR) THEN BEGIN
        plot,vc_x,vc_y[0,*],psym=10,xtitle="Radius [kpc]",ytitle="Velocity [km/s]",yrange = [0,maxvc],thick = thicks[0];yrange = [0,260]
        oplot,vc_x,vc_y[0,*],psym=10,color=colors[0],thick = thicks[0]
    ENDIF ELSE plot,vc_x,vc_y[0,*],psym=10,color=colors[0],xtitle="Radius [kpc]",ytitle="Velocity [km/s]",yrange = [0,60];yrange = [0,260]
ENDIF ELSE BEGIN 
;        window,9
    plot,vc_x,vc_y[0,*],psym=10,xtitle="Radius [kpc]",ytitle="Velocity [km/s]",yrange = [0,maxvc],thick = thicks[0] ;yrange = [0,260];,title='"Line Profiles" -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years'
ENDELSE
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    oplot,vc_x,vc_y[iFile,*],linestyle=linestyles[0],psym=10,color=colors[iFile],thick = thicks[iFile]
ENDFOR
IF (N_ELEMENTS(files) gt 1) THEN legend,keys,color = colors,thick = thicks,/right,linestyle = fltarr([N_ELEMENTS(files)]),/bottom,charsize = l_charsize
;IF (KEYWORD_SET(outfile)) THEN device,/close
;ENDIF

;********************************* Line Profile ***********************************************************
;IF (show_line) THEN BEGIN
maxl = 5e7 ;1e8
IF (KEYWORD_SET(outfile)) THEN BEGIN
;   device,filename=outfile+'line_.eps',/color,bits_per_pixel= 8,/times 
    IF KEYWORD_SET(COLOR) THEN BEGIN
        plot,v_bin,v_y[1,0,*],psym=10,xtitle="Velocity",ytitle="Gas Mass",yrange = [0,maxl],thick = thicks[0]
    ENDIF ELSE plot,v_bin,v_y[1,0,*],psym=10,color=colors[0],xtitle="Velocity",ytitle="Gas Mass",yrange = [0,maxl],thick = thicks[0]
ENDIF ELSE BEGIN 
;   window,4
    plot,v_bin,v_y[1,0,*],psym=10,xtitle="Velocity",ytitle="Gas Mass",yrange = [0,maxl],thick = thicks[0];title='"Line Profiles" -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years'
ENDELSE
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    oplot,v_bin,v_y[1,iFile,*],linestyle=linestyles[0],psym=10,color=colors[iFile],thick = thicks[iFile]
    oplot,v_bin,v_y[0,iFile,*],linestyle=linestyles[1],psym=10,color=colors[iFile],thick = thicks[iFile]
ENDFOR
legend,["HI "+legend_t,"H2 "+legend_t],linestyle=linestyles,charsize = l_charsize;,color=colors,thick = thicks
;IF (N_ELEMENTS(files) gt 1) THEN legend,keys,color = colors,thick = thicks,/right,charsize = l_charsize,linestyle = fltarr([N_ELEMENTS(files)])
;IF (KEYWORD_SET(outfile)) THEN device,/close
;ENDIF

;********************************* Temperature Distro ***********************************************************
;IF (show_temp) THEN BEGIN
IF (KEYWORD_SET(outfile)) THEN BEGIN
;    device,filename=outfile+'temp_.eps',/color,bits_per_pixel= 8,/times 
    IF KEYWORD_SET(COLOR) THEN BEGIN
        plot,t_bin,t_y[1,0,*],psym=10,thick = thicks[0],xtitle="Temperature",ytitle="Gas Mass",yrange = [1e3,1e8],/ylog
;           oplot,t_bin,t_y[1,0,*],psym=10,color=colors[0],thick = thicks[0]
    ENDIF ELSE plot,t_bin,t_y[1,0,*],psym=10,color=colors[0],thick = thicks[0],xtitle="Temperature",ytitle="Gas Mass",yrange = [1e3,1e8],/ylog
ENDIF ELSE BEGIN 
;        window,6      
    plot,t_bin,t_y[1,0,*],psym=10,xtitle="Temperature",ytitle="Gas Mass",yrange = [1e3,1e8],thick = thicks[0],/ylog ;,title='Temperature Distribution -- '  + strtrim(i*dt*dDelta*timeunit,2) + ' years';,yrange=[0,20000]
ENDELSE
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    oplot,t_bin,t_y[1,iFile,*],psym=10,linestyle=linestyles[0],color=colors[iFile],thick = thicks[iFile]
    oplot,t_bin,t_y[0,iFile,*],psym=10,linestyle=linestyles[1],color=colors[iFile],thick = thicks[iFile]
ENDFOR
;legend,["HI "+legend_t,"H2 "+legend_t],linestyle=linestyles,charsize = l_charsize ;,color=colors,thick = thicks
;IF (N_ELEMENTS(files) gt 1) THEN legend,keys,color = colors,thick = thicks,/bottom,charsize = l_charsize,linestyle = fltarr([N_ELEMENTS(files)])
;IF (KEYWORD_SET(outfile)) THEN device,/close        
;ENDIF
 
;********************************* Density Distro ***********************************************************
;IF (show_den) THEN BEGIN    
IF (KEYWORD_SET(outfile)) THEN BEGIN
;    device,filename=outfile+'dens_.eps',/color,bits_per_pixel= 8,/times 
    IF KEYWORD_SET(COLOR) THEN BEGIN
        plot,rho_bin,rho_y[1,0,*],psym=10,thick = thicks[0],xtitle="Density",ytitle="Gas Mass" ,yrange = [1e3,1e8],/ylog
;           oplot,rho_bin,rho_y[1,0,*],psym=10,color=colors[0],thick = thicks[0]
    ENDIF ELSE plot,rho_bin,rho_y[1,0,*],psym=10,color=colors[0],thick = thicks[0],xtitle="Density",ytitle="Gas Mass" ,yrange = [1e3,1e8],/ylog
ENDIF ELSE BEGIN 
;        window,7
    plot,rho_bin,rho_y[1,0,*],psym=10,xtitle="Density",ytitle="Gas Mass" ,yrange = [1e3,1e8],thick = thicks[0],/ylog ;,title='Density Distribution -- '  + strtrim(i*dt*dDelta*timeunit,2) + ' years'
ENDELSE
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    oplot,rho_bin,rho_y[1,iFile,*],psym=10,linestyle=linestyles[0],color=colors[iFile],thick = thicks[iFile]
    oplot,rho_bin,rho_y[0,iFile,*],psym=10,linestyle=linestyles[1],color=colors[iFile],thick = thicks[iFile]
ENDFOR
;legend,["HI "+legend_t,"H2 "+legend_t],linestyle=linestyles,charsize= l_charsize ;,color=colors,thick = thicks
;IF (N_ELEMENTS(files) gt 1) THEN legend,keys,color = colors,thick = thicks,/right,charsize = l_charsize,linestyle = fltarr([N_ELEMENTS(files)])
;IF (KEYWORD_SET(outfile)) THEN device,/close
;ENDIF

;********************************* Radial Profile ***********************************************************
;IF (show_prof) THEN BEGIN
maxprof = 7
IF (KEYWORD_SET(outfile)) THEN BEGIN
;    device,filename=outfile+'profile_.eps',/color,bits_per_pixel=
;    8,/times 
    IF KEYWORD_SET(COLOR) THEN BEGIN
        plot,r_bin,prof_y[1,0,*],psym=10,xtitle="Radius [kpc]",ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,/ylog,thick = thicks[0],xrange =[0,maxprof],yrange =[1,200]
;           oplot,r_bin,prof_y[1,0,*],psym=10,color=colors[0],thick = thicks[0]
    ENDIF ELSE plot,r_bin,prof_y[1,0,*],psym=10,color=colors[0],xtitle="Radius [kpc]",ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,/ylog,thick = thicks[0],xrange =[0,maxprof],yrange =[1,200]
ENDIF ELSE BEGIN 
;        window,5
    plot,r_bin,prof_y[1,0,*],psym=10,xtitle="Radius [kpc]",ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,/ylog,thick = thicks[0],xrange =[0,maxprof],yrange =[1,200] ;,title='Gas Surface Density Profile -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years',/ylog,yrange = [1e2,3e7]
ENDELSE
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    oplot,r_bin,prof_y[1,iFile,*],psym=10,linestyle=linestyles[0],color=colors[iFile],thick = thicks[iFile]
    oplot,r_bin,prof_y[0,iFile,*],psym=10,linestyle=linestyles[1],color=colors[iFile],thick = thicks[iFile]
    oplot,prof_stellar[0,iFile,*],prof_stellar[1,iFile,*]/1e6,psym=10,color=colors[iFile],thick = thicks[iFile],linestyle = 3
ENDFOR
;legend,["HI "+legend_t,"H2 "+legend_t],linestyle=linestyles,/right,charsize=l_charsize ;,color=colors,thick = thicks
;;;;;legend,STRMID(STRTRIM(halo_stats[*].HImass + halo_stats[*].H2mass,2),0,4) +STRMID(STRTRIM(halo_stats[*].HImass + halo_stats[*].H2mass,2),3,/REVERSE_OFFSET) + ", " +  STRMID(STRTRIM(halo_stats[*].smass,2),0,4) +STRMID(STRTRIM(halo_stats[*].smass,2),3,/REVERSE_OFFSET),color = colors,thick = thicks,/top,/right,linestyle = fltarr([N_ELEMENTS(files)]),charsize = l_charsize
;IF (KEYWORD_SET(outfile)) THEN device,/close
;ENDIF

;******************************* Stellar Profile ****************************************

;IF (stellar_prof) THEN BEGIN
;IF (KEYWORD_SET(outfile)) THEN BEGIN
;    device,filename=outfile+'stellarprofile_.eps',/color,bits_per_pixel= 8,/times 
;    IF KEYWORD_SET(COLOR) THEN BEGIN
;        plot,prof_stellar[0,1,*],prof_stellar[1,1,*],psym=10,xtitle="Radius [kpc]",yrange = [1e3,1e9],ytitle="Stellar Mass Surface Density [M"+sunsymbol()+textoidl(' kpc^2')+"]" ,/ylog,thick = thicks[0]
;           oplot,r_bin,prof_y[1,0,*],psym=10,color=colors[0],thick = thicks[0]
;    ENDIF ELSE plot,prof_stellar[0,1,*],prof_stellar[1,1,*],psym=10,color=colors[0],xtitle="Radius [kpc]",yrange = [1e3,1e9],ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' kpc^2')+"]" ,/ylog 
;ENDIF ELSE BEGIN 
;        window,9
;    plot,prof_stellar[0,0,*],prof_stellar[1,0,*],psym=10,xtitle="Radius [kpc]",ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' kpc^2')+"]" ,/ylog,thick = thicks[0],yrange = [1e3,1e9] ;,title='Gas Surface Density Profile -- ' + strtrim(i*dt*dDelta*timeunit,2) + ' years',/ylog,yrange = [1e2,3e7]
;ENDELSE
;FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
;    oplot,prof_stellar[0,iFile,*],prof_stellar[1,iFile,*],psym=10,color=colors[iFile],thick = thicks[iFile]
;ENDFOR
legend,["HI "+legend_t,"H2 "+legend_t],linestyle=linestyles,/right,charsize=l_charsize,/bottom
;;,color=colors,thick = thicks
legend,STRMID(STRTRIM(halo_stats[*].smass,2),0,4) +STRMID(STRTRIM(halo_stats[*].smass,2),3,/REVERSE_OFFSET) ,color = colors,thick = thicks,/top,/right,linestyle = fltarr([N_ELEMENTS(files)]),charsize = l_charsize
;IF (KEYWORD_SET(outfile)) THEN device,/close
;ENDIF

;*************************** Schmidt-Kennicutt **************************************
xsigmalow = findgen(300)/100 - 1.
xsigma = 10.0^(findgen(600)/100 - 3.) ;xsigmalow
ysigma=2.5e-4*xsigma^1.4
ysigma1=2.5e-4*xsigma^(1.4 + 0.15)
ysigma2=2.5e-4*xsigma^(1.4 - 0.15)
ysigmalow = xsigmalow*2.4 - 5.0

IF (KEYWORD_SET(outfile)) THEN BEGIN
;    device,filename=outfile+'sk.eps',/color,bits_per_pixel= 8,/times,xsize = 16,ysize = 12
;    !X.charsize = 0.75
    plot,alog10(xsigma),alog10(ysigma),xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,  xrange = [-3,3], yrange = [-6,3]
ENDIF ELSE BEGIN
;    window,5,xsize = 712,ysize = 392
    plot,alog10(xsigma),alog10(ysigma),xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,  xrange = [-3,3], yrange = [-6,3]
ENDELSE
;oplot,xsigmalow,ysigmalow
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    IF KEYWORD_SET(COLOR) THEN $
      oplot,[alog10(halo_stats[iFile].gasSD),alog10(halo_stats[iFile].gasSD)],[alog10(halo_stats[iFile].starSD),alog10(halo_stats[iFile].starSD)],psym=4,color=colors[iFile],thick = thicks[iFile] ELSE BEGIN
        oplot,[alog10(halo_stats[iFile].gasSD),alog10(halo_stats[iFile].gasSD)],[alog10(halo_stats[iFile].starSD),alog10(halo_stats[iFile].starSD)],psym=symbols[ifile],color=colors[iFile],thick = thicks[iFile]
        legend,keys,color=colors,thick = thicks,charsize = l_charsize,linestyle = fltarr([N_ELEMENTS(files)]),psym = symbols,/right
    ENDELSE
ENDFOR
;;;;legend,keys,color=colors,thick = thicks,charsize = l_charsize,linestyle = fltarr([N_ELEMENTS(files)]),/right
IF (KEYWORD_SET(outfile)) THEN BEGIN
    device,/close
;    !X.charsize = 1        
ENDIF

!p.multi = 0
IF (KEYWORD_SET(outfile)) THEN BEGIN
    device,filename=outfile+'sk.eps',/color,bits_per_pixel= 8,/times,xsize = 16,ysize = 12
;    !X.charsize = 0.75
    plot,alog10(xsigma),alog10(ysigma),xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,  xrange = [-3,3], yrange = [-6,3]
ENDIF ELSE BEGIN
    window,5,xsize = 712,ysize = 392
    plot,alog10(xsigma),alog10(ysigma),xtitle=textoidl('Log \Sigma')+"!lgas!n [M"+sunsymbol()+" pc!u-2!n]",ytitle=textoidl('Log \Sigma')+"!lSFR!n [M"+sunsymbol()+" kpc!u-2!n yr!u-1!n]", xstyle=1, ystyle=1,  xrange = [-3,3], yrange = [-6,3]
ENDELSE
;oplot,xsigmalow,ysigmalow
FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
    IF KEYWORD_SET(COLOR) THEN $
      oplot,[alog10(halo_stats[iFile].gasSD),alog10(halo_stats[iFile].gasSD)],[alog10(halo_stats[iFile].starSD),alog10(halo_stats[iFile].starSD)],psym=4,color=colors[iFile],thick = thicks[iFile] ELSE BEGIN
      oplot,[alog10(halo_stats[iFile].gasSD),alog10(halo_stats[iFile].gasSD)],[alog10(halo_stats[iFile].starSD),alog10(halo_stats[iFile].starSD)],psym=symbols[iFile],color=colors[iFile],thick = thicks[iFile]      
  ENDELSE
ENDFOR
legend,keys,color=colors,thick = thicks,charsize = l_charsize,linestyle = fltarr([N_ELEMENTS(files)]),psym = symbols
IF (KEYWORD_SET(outfile)) THEN BEGIN
    device,/close
;    !X.charsize = 1        
ENDIF



!p.multi = [0,2,1]
;***************** Baryonic Tully-Fisher
IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+'_btf.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2 ELSE window,5,xsize = 712,ysize = 392
;legend,STRMID(STRTRIM(halo_stats[*].HImass + halo_stats[*].H2mass,2),0,4) +STRMID(STRTRIM(halo_stats[*].HImass + halo_stats[*].H2mass,2),3,/REVERSE_OFFSET) + ", " +  STRMID(STRTRIM(halo_stats[*].smass,2),0,4) +STRMID(STRTRIM(halo_stats[*].smass,2)

readcol,'/astro/users/christensen/code/MolecH/BTF.dat',name,V_f,M_star,M_gas,mu_0,R_d,B_V,gamma_max,gamma_pop,gamma_acc,FORMAT = '(A9D)'
M_star = M_star*1e10
M_gas = M_gas*1e10
M_disk = M_star + M_gas
plot,M_disk,M_star/M_disk,psym = 1,xtitle = 'Disk Mass [Msol]',ytitle = 'Stellar Fraction',/xlog,xrange= [1e8,1e12],yrange = [-0.2,1.2]
IF KEYWORD_SET(COLOR) THEN BEGIN
    FOR i = 0, N_ELEMENTS(halo_stats) - 1 DO oplot,[halo_stats[i].HImass + halo_stats[i].smass,halo_stats[i].HImass + halo_stats[i].smass],[halo_stats[i].smass/(halo_stats[i].HImass + halo_stats[i].smass),halo_stats[i].smass/(halo_stats[i].HImass + halo_stats[i].smass)],color = colors[i],psym = 4,SYMSIZE = 1.5
    legend, keys, color = colors, thick = thicks, /top, /left, linestyle = fltarr([N_ELEMENTS(files)]),charsize = l_charsize
ENDIF ELSE BEGIN
    FOR i = 0, N_ELEMENTS(halo_stats) - 1 DO oplot,[halo_stats[i].HImass + halo_stats[i].smass,halo_stats[i].HImass + halo_stats[i].smass],[halo_stats[i].smass/(halo_stats[i].HImass + halo_stats[i].smass),halo_stats[i].smass/(halo_stats[i].HImass + halo_stats[i].smass)],color = colors[i],psym = symbols[i],SYMSIZE = 1.5
    legend, keys, color = colors, thick = thicks, /top, /left, linestyle = fltarr([N_ELEMENTS(files)]),psym = symbols, charsize = l_charsize
ENDELSE

plot,V_f,M_star + M_gas,psym = 1,/xlog,/ylog,xrange = [40,500],yrange = [1e8,4e11],xtitle = textoidl('V_{200} [km/s]'),ytitle = textoidl('Baryonic Mass [M')+sunsymbol()+']'
IF KEYWORD_SET(COLOR) THEN $
  FOR i = 0, N_ELEMENTS(halo_stats) - 1 DO oplot,[halo_stats[i].v_f,halo_stats[i].v_f], [halo_stats[i].smass + halo_stats[i].HImass,halo_stats[i].smass + halo_stats[i].HImass],color = colors[i],psym = 4,SYMSIZE = 1.5 ELSE $
  FOR i = 0, N_ELEMENTS(halo_stats) - 1 DO oplot,[halo_stats[i].v_f,halo_stats[i].v_f], [halo_stats[i].smass + halo_stats[i].HImass,halo_stats[i].smass + halo_stats[i].HImass],color = colors[i],psym = symbols[i],SYMSIZE = 1.5

!p.multi = 0
;********************************* SFR Diagram ***********************************************************
IF (show_sfr) THEN BEGIN
    maxsfr = 0.4;10 ;0.4;MAX(sfr_y[0,*,*])
    IF (KEYWORD_SET(outfile)) THEN BEGIN
        device,filename=outfile+'_sfr.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset =  2 
;        !X.charsize = 0.75
        plot,sfr_bin,sfr_y[0,0,*],psym=10,xtitle="Time [Gyr]",ytitle='SFR [M'+sunsymbol()+'/yr]',yrange = [0,maxsfr],thick = thicks[0];,charsize = 2.0;thicks[0]
    ENDIF ELSE BEGIN
        window,4,xsize = 712,ysize = 392
        plot,sfr_bin,sfr_y[0,0,*],psym=10,xtitle="Time [Gyr]",ytitle='SFR [M'+sunsymbol()+'/yr]',yrange = [0,maxsfr],thick = thicks[0] ;,title='Density Distribution -- '  + strtrim(i*dt*dDelta*timeunit,2) + ' years'
    ENDELSE
    FOR ifile = 0, N_ELEMENTS(files) - 1 DO oplot,sfr_bin,sfr_y[0,iFile,*],psym=10,linestyle=linestyles[0],color=colors[iFile],thick = thicks[iFile]
    legend,keys,color=colors,thick = thicks,charsize = l_charsize,linestyle = fltarr([N_ELEMENTS(files)])
    IF (KEYWORD_SET(outfile)) THEN BEGIN
        device,/close
;        !X.charsize = 1        
    ENDIF
ENDIF

;!p.multi = 0
;IF (KEYWORD_SET(outfile)) THEN device,filename=outfile+'_profs.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 9,xoffset =  2,yoffset =  2
; plot,r_bin,prof_y[1,0,*],psym=10,xtitle="Radius [kpc]",yrange = [1,1000],ytitle="Gas Mass Surface Density [M"+sunsymbol()+textoidl(' pc^2')+"]" ,/ylog,thick = thicks[0],xrange = [0,4]
;FOR ifile = 0, N_ELEMENTS(files) - 1 DO BEGIN
;    oplot,r_bin,prof_y[1,iFile,*],psym=10,linestyle=linestyles[0],color=colors[iFile],thick = thicks[iFile]
;    oplot,r_bin,prof_y[0,iFile,*],psym=10,linestyle=linestyles[1],color=colors[iFile],thick = thicks[iFile]
;    oplot,prof_stellar[0,iFile,*],prof_stellar[1,iFile,*]/1e6,psym=10,color=colors[iFile],thick = thicks[iFile],linestyle = 1
;ENDFOR
;legend,["HI","H2", "Stars"],linestyle=[0,2,1],/right,charsize=1.2,/bottom
;;;,color=colors,thick = thicks
;legend,keys ,color = colors,thick = thicks,/top,/right,linestyle = fltarr([N_ELEMENTS(files)])
;IF (KEYWORD_SET(outfile)) THEN device,/close 
IF NOT keyword_set(outfile) THEN stop
END
 
