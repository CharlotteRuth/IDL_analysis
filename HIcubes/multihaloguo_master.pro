PRO LGanalysis
;stepst = ['00504','00512']
;outplot = '/astro/net/nbody1/christensen/Schmidtlaw/data/LG.504.512.guo.ps'
;psym = [4,5]
;color = [50,240]

stepst = ['00512']
outplot = '/astro/net/nbody1/christensen/Schmidtlaw/data/LG.512.guo.ps'
psym = [4]
color = [240]

files = '/astro/net/scratch2/christensen/LG/gas/steps/LG.'+stepst+'.dir/LG.'+stepst


stepst = ['00512']
outplot = '/astro/net/nbody1/christensen/Schmidtlaw/data/LG.512.guo.ps'
psym = [4]
color = [240]

files = ['/astro/net/scratch2/christensen/LG/gas/steps/LG.00480.dir/LG.00480','/astro/net/scratch2/christensen/LG/gas_oldfeedback/steps/LG.00480.dir/LG.00480']
outplot = '/astro/net/nbody1/christensen/Schmidtlaw/data/LG.480.fb.guo.ps'
psym = [4,4]
color = [240,50]
multiHalGuo, files, color = color, psym = psym,imf = 2,key = ['LG.00480','LG.00480 old fb'],outplot = outplot
;multiHaloGuo, files, color = color, psym = psym,imf = 2,key = stepst,/kcorrect,outplot = outplot
END

PRO h516analysis
dir = '/astro/net/'
;step = '00168'
;files = ['h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.'+step+'.dir/h516.cosmo25cmb.1536g1MBWK.'+step, $
;         'h516.cosmo25cmb.1536g2HBWK_nocool/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step, $
;         'h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step]
step = '00512'
files = ['scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.'+step+'.dir/h516.cosmo25cmb.1536g1MBWK.'+step, $
         'scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step,$
          'scratch1/cbrook/sims/h516.cosmo25cmb.2304g4bwdK/h516.cosmo25cmb.2304g4bwdK.'+step+'/h516.cosmo25cmb.2304g4bwdK.'+step]
key = ['No H2', 'H2','No H2, 2304']
color =  [50,240,70]
psym = [4,4,4]
outplot = '/astro/net/nbody1/christensen/Schmidtlaw/data/h516.H2.noH2.2304.guo.ps'
;step = '00168'
;files = ['h516.cosmo25cmb.1536g2HBWK_nocool/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step, $
;         'h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.'+step+'.dir/h516.cosmo25cmb.1536g2HBWK.'+step]
;key = ['H2, w/o cooling', 'H2 w/ cooling']
 
step = '00384'
dir =  "/astro/net/scratch2/christensen/MolecH/Cosmo/"
files = ['h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00384.dir/h516.cosmo25cmb.1536g3HBWK.00384','h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00384.dir/h516.cosmo25cmb.1536g6MbwK.00384']
key = ['H2','no H2']
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/guo.eps'
color = [120,240]
psym = [4,4]


dir = '/astro/net/scratch2/christensen/MolecH/Cosmo/'
files = ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.00512/h516.cosmo25cmb.3072g1MBWK.00512', $
         'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/steps/h516.cosmo25cmb.1536g6MbwK.00512.dir/h516.cosmo25cmb.1536g6MbwK.00512', $
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps_noH2SF/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512.dir/h516.cosmo25cmb.1536g3HBWK_noH2SF.00512', $
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir/h516.cosmo25cmb.1536g3HBWK.00512', $
        'h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/steps/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464.dir/h516.cosmo25cmb.1536g6HBWK.jeans.prev.00464']
color = [0,30,90,140,240]
psym = [4,4,4,4,4]
key = ['High Res','no H2','H2','H2 SF','MC SF']
symsize = [1,1,1,1,1]
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.guo.eps'
files = dir + files

files = [ '/astro/net/nbody1/abrooks/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.00324/h516.cosmo25cmb.3072g1MBWK.00324',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00320.dir/h516.cosmo25cmb.3072g14HBWK.00320']
key = ['h516, no H2','h516, H2']
outplot = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/Twins.3072g.guo.324.kcorrect.eps']
color = [120,240]
psym = [4,4]
haloids = [[1,2,4,6,15,28],[1,2,4,6,15,28]]
outplot = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.guo.324.eps'

multiHaloGuo,files,color = color,psym = psym,key = key,symsize=symsize,xrange = [5,12],yrange = [4,10.5],haloids = haloids,halomasses = halomasses,/kcorrect,outplot = outplot;,xrange = [5,11],yrange = [4,9.5]

outplot = 0
IF KEYWORD_SET(outplot) THEN BEGIN
    device, filename='/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g.guoratio.324.kcorrect.eps',/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
ENDIF 

plot,alog10(halomasses[*,0]),alog10(halomasses[*,1]),psym = 5,xtitle = 'Log(Mass Halo)',ytitle = 'Log(Mstar_H2/Mstar_{no H2})',yrange = [-0.5,3],xrange = [8.4,10.6]
stop
IF KEYWORD_SET(outplot) THEN device,/close
END

PRO h516analysis_z
prefix = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/'
base = 'h516.cosmo25cmb.3072g1MBWK'
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK.guo.eps'
bad = [168,204,264]
start = 36

;prefix = '/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/'
;base = 'h603.cosmo50cmb.3072gs1MbwK'
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g.guo.eps'
;bad = [204,312,372,384,504]
;;start = 120
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g.guo.kcorrect.eps'
;bad = [204,276,312,372,384,504]
;start = 132;120

imf = 1

last = 512; 168;512
dt = 12
nsteps = last/dt

files = strarr(nsteps - N_ELEMENTS(bad) - start/dt + 1 + 1)
key = strarr(nsteps - N_ELEMENTS(bad) - start/dt + 1 + 1)
color = fltarr(nsteps - N_ELEMENTS(bad) - start/dt + 1 + 1)
psym = fltarr(nsteps - N_ELEMENTS(bad) - start/dt + 1 + 1) + 1
symsize = fltarr(nsteps - N_ELEMENTS(bad) - start/dt + 1 + 1) + 1
steps = strarr(nsteps - N_ELEMENTS(bad) - start/dt + 1 + 1)
it = 0

FOR i = start/dt, nsteps DO BEGIN 
    step = i*dt
    if where(step eq bad) eq -1  THEN BEGIN
        if (step lt 10) THEN step = '0000'+STRTRIM(step,2) ELSE BEGIN
            if (step lt 100) THEN step = '000'+STRTRIM(step,2) ELSE step = '00'+STRTRIM(step,2)
        ENDELSE
        steps[it] = step
;        filename = [prefix + base + '.' + step + '/' + base + '.' + step]
;        key = [step]
;        multiHaloGuo,[filename],key = key,color = 240,yrange = [4,11],xrange = [7,13] ;,outplot = outplot
        files[it] = prefix + base + '.' + step + '/' + base + '.' + step
        key[it] = step
        color[it] = 256/nsteps*i
        it = it + 1
    ENDIF
ENDFOR
step = '00512'
steps[it] = step
files[it] = prefix + base + '.' + step + '/' + base + '.' + step
key[it] = step
color[it] = 1
psym[it] = 2
symsize[it] = 1.5
;multiHaloGuo,files,color = color,yrange = [4,11],xrange = [7,13],psym = psym,symsize = symsize,imf = imf,outplot = outplot;,/kcorrect;,key = key,

FOR i = 0, N_ELEMENTS(files) - 1  DO BEGIN 
     multiHaloGuo,files[0:i],color = color[0:i],yrange = [4,11],xrange = [7,13],psym = psym[0:i],symsize = symsize[0:i],imf = imf,outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK.guo_movie/h516.cosmo25cmb.3072g1MBWK.'+steps[i]+'.guo.eps',/kcorrect ;,key = key
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g.guo_movie/h603.cosmo50cmb.3072gs1MbwK.'+steps[i]+'.guo.eps'
 ENDFOR

END

PRO h516analysis_zsteps
prefix = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'

;base = 'h516.cosmo25cmb.1536g3HBWK'
;outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK.guo.eps'
;steps = ['00015','00030','00045','00060','00075','00090','00105','00120','00135','00150','00165','00180','00195','00210','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00468','00480','00492','00504','00512']

base = 'h516.cosmo25cmb.1536g6MbwK'
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6MbwK/h516.cosmo25cmb.1536g6MbwK.guo.eps'
steps = ['00192','00240','00328','00384','00406','00455','00480','00512']

base = 'h516.cosmo25cmb.1536g6HBWK'
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g6HBWK/Jeans_oldLW/h516.cosmo25cmb.1536g6HBWK.guo.eps'
steps = ['00012','00024','00036','00048','00060','00072','00084','00096','00108','00132','00144','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324','00336','00348','00360','00372','00384','00396','00408','00420','00432','00444','00456','00464'] ;120,168

imf = 1
nsteps = N_ELEMENTS(steps);last/dt
files = strarr(nsteps)
key = strarr(nsteps)
color = fltarr(nsteps)
psym = fltarr(nsteps) + 1
symsize = fltarr(nsteps) + 1
it = 0

FOR i = 0, nsteps - 1 DO BEGIN 
    step = steps[i]
    ;if where(step eq bad) eq -1  THEN BEGIN
;        filename = [prefix + base + '.' + step + '/' + base + '.' + step]
;        key = [step]
;        multiHaloGuo,[filename],key = key,color = 240,yrange = [4,11],xrange = [7,13] ;,outplot = outplot
        files[it] = prefix + base + '/steps/' + base + '.' + step + '.dir/' + base + '.' + step
        files[it] = prefix + base + '/Jeans_oldLW/steps/' + base + '.jeans.prev.' + step + '.dir/' + base + '.jeans.prev.' + step
        key[it] = step
        color[it] = 256.0/MAX(fix(steps))*fix(step)
        it = it + 1
 ;   ENDIF
ENDFOR
color[it - 1] = 1
psym[it - 1] = 2
symsize[it - 1] = 1.5
multiHaloGuo,files,color = color,yrange = [4,11],xrange = [7,13],psym = psym,symsize = symsize,imf = imf,outplot = outplot;,/kcorrect;,key = key,
END



PRO h603analysis
files = ['/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g3HBWK/steps/h603.cosmo50cmb.1536g3HBWK.00144.dir/h603.cosmo50cmb.1536g3BWK.00144', $
        '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK.BUG/h603.cosmo50cmb.2304g2bwK.00144/h603.cosmo50cmb.2304g2bwK.00144', $
        '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG/h603.cosmo50cmb.2304g5bwK.00144.dir/h603.cosmo50cmb.2304g5bwK.00144']
key = ['H2','low threshhold','high threshold']
outplot = '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g3HBWK/steps/h603.cosmo50cmb.1536g3HBWK.00144.dir/h603.cosmo50cmb.1536g3BWK.00144.guo.ps'
color = [140,240,50]
psym = [4,4,4]

files = ['/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.00072/h603.cosmo50cmb.3072gs1MbwK.00072',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00072.dir/h603.cosmo50cmb.3072g14HBWK.00072']
key = ['no H2','H2']
outplot = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g.guo.072.eps']

files = ['/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.00324/h603.cosmo50cmb.3072gs1MbwK.00324',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324']
key = ['no H2','H2']
outplot = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g.guo.324.kcorrect.eps']

;files = [ '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324']
;key = ['H2']
;outplot = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK.guo.324.eps']
haloids = [[1,2,3,4,5,6,15,16,18,20],$
           [1,2,3,4,5,7,14,16,18,19]]
multiHaloGuo,files,/color,key = key,xrange = [5,12],yrange = [4,10.5],haloids = haloids,halomasses = halomasses,/kcorrect;,halomasses = halomasses
IF KEYWORD_SET(outplot) THEN BEGIN
    device, filename='/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g.guoratio.324.kcorrect.eps',/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
ENDIF 
plot,alog10(halomasses[*,0]),alog10(halomasses[*,1]),psym = 5,xtitle = 'Log(Mass Halo)',ytitle = 'Log(Mstar_H2/Mstar_{no H2})',yrange = [-1,1],xrange = [9,11.5]
IF KEYWORD_SET(outplot) THEN device,/close
stop
END

;outplot ='/astro/users/christensen/plots/h603/h603.cosmo50cmb.3072g14HBWK.zsteps_guo.eps'
;outplot ='/astro/users/christensen/plots/h603/h603.cosmo50cmb.3072g14HBWK.zsteps_guo.fabio.eps'
PRO h603analysis_zsteps,outplot = outplot

IF NOT KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'x'
    formatplot
    window,0
ENDIF

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/'
base = 'h603.cosmo50cmb.3072g14HBWK'
steps =   ['00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324']
z =       ['3.436','2.999','2.654','2.372','2.138','1.939','1.766','1.616','1.483','1.364','1.257','1.160','1.072','0.991','0.916','0.846','0.782','0.722','0.665','0.612','0.562','0.515']
tags =    ['3.5'  ,    ' ',    ' ',    ' ',    ' ',    '2',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    '1',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',  '0.5']
haloids = [[    2],[    2],[    2],[    2],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1]]
imf = 1
nsteps = N_ELEMENTS(steps);last/dt
files = strarr(nsteps)
key = strarr(nsteps)
color = fltarr(nsteps)
psym = fltarr(nsteps) - 1
symsize = fltarr(nsteps) + 1
xrange = [8,12]
yrange = [5,11]

FOR i = 0, nsteps - 1 DO BEGIN 
    step = steps[i]
    files[i] = prefix + base + '/steps/' + base + '.' + step + '.dir/' + base + '.' + step
    key[i] = step
ENDFOR
multiHaloGuo,files,yrange = yrange,xrange = xrange,psym = psym,imf = imf,haloids = haloids,tags = tags;,/kcorrect;,key = key,
stop
END

PRO multianalysis,outplot = outplot2,color = color, kcorrect = kcorrect,moster = moster
;outplot = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/Twins.3072g.guo.324.kcorrect.eps']
;outplot = '/astro/users/christensen/plots/h603/stellarmass_halomass'
;outplot = '/astro/users/christensen/plots/h603/stellarmass_halomass.fabio'

;at least 40,000 DM particles
formatplot,outplot = outplot2
IF KEYWORD_SET(outplot2) THEN BEGIN
    fgcolor = 0
    bgcolor = 255
    if keyword_set(kcorrect) THEN outplot = outplot2 + 'obs' else outplot = outplot2
    device, filename=outplot + '_guoz0.eps',/COLOR,bits_per_pixel= 8,/times,ysize=12,xsize=18;,/inch
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    window,0
ENDELSE
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
;USERSYM, COS(A), SIN(A), /FILL
filesh603 = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.00512.dir/h603.cosmo50cmb.3072g1bwK.00512',$
;             '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.3072g1bwK/steps/h603.cosmo50cmb.3072g1bwK.00512.dir/h603.cosmo50cmb.3072g1bwK.00512',$
;             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.00512.dir/h603.cosmo50cmb.2304g2bwK.00512',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.00512.dir/h603.cosmo50cmb.3072gs1MbwK.00512',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00512.dir/h603.cosmo50cmb.3072g14HBWK.00512']
haloidsh603 = [[1,2,3,8,12,1],$
               [1,2,3,7,8,16],$
               [1,2,3,8,7,13]]

filesh986 = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072g1bwK/h986.cosmo50cmb.3072g1bwK.00512/h986.cosmo50cmb.3072g1bwK.00512',$
               '/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.00512/h986.cosmo50cmb.3072gs1MbwK.00512',$
               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir/h986.cosmo50cmb.3072g14HBWK.00512']
haloidsh986 = [[1, 2, 4, 9, 11, 13, 14],$
               [1, 2, 3, 4,  5,  7, 10],$
               [1, 2, 3, 9, 11, 14, 17]]

;filesh516 = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492',$
;             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512']
;haloidsh516 = [[1,2,3,4,6,13],$
;               [1,2,3,4,6,13]]

filesh516 = [ '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g3BWK/steps/h516.cosmo25cmb.3072g3BWK.00512.dir/h516.cosmo25cmb.3072g3BWK.00512',$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492' ,$
                  '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512']
haloidsh516 = [[1,2,4,6,16],$
               [1,2,3,6,13],$
               [1,2,3,6,13]]     

filesh799 = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3BWK/steps/h799.cosmo25cmb.3072g3BWK.00512.dir/h799.cosmo25cmb.3072g3BWK.00512',$
;'/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g3bwdK/steps/h799.cosmo25cmb.3072g3bwdK.00512.dir/h799.cosmo25cmb.3072g3bwdK.00512',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/steps/h799.cosmo25cmb.3072g1MBWK.00512.dir/h799.cosmo25cmb.3072g1MBWK.00512', $
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00512.dir/h799.cosmo25cmb.3072g14HBWK.00512']
haloidsh799 = [[1,6,11,15,16],$
               [1,6,11,15,16],$
               [1,6,10,15,16]]


;xrange = [8.5,12]
xrange = [9,11.5]
;yrange = [5,11]
yrange = [5,10.5]

;xrange = [8.4,10.65]
;yrange = [0,13]

n = N_ELEMENTS(filesh603)
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    obscolor = 0 ;fgcolor
    if not keyword_set(colors) then  colors  = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(psym) THEN psym = fltarr(n) + 4

    obscolor = 50
    obssym = 15
;    fgcolor = 255
    colors = [fgcolor,fgcolor,245]
    psym = [4,5,16]

;fgcolor = 0
    obscolor = 50
    obssym = 2
    colors = [120,fgcolor,254]
    psym = [4,6,16]
    ctables = [0,0,39]
ENDIF ELSE BEGIN
    loadct,0    
    obscolor = fgcolor
    obssym = 15
    colors = [100,100,fgcolor];(fltarr(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    colors = [0,0,fgcolor]
    ctables = [0,0,0]
    IF NOT KEYWORD_SET(psym) THEN  psym = (findgen(n)+2)*2
ENDELSE
psymh603  = psym      ;[4,4,4]
colorh603 = colors     ;[120,240,50]
psymh986  = psym      ;[4,4,4]
colorh986 = colors     ;[120,240,50]
psymh516  = psym;[1:2] ;[5,5]
colorh516 = colors;[1:2];[240,50]
psymh799  = psym      ;[4,4,4]
colorh799 = colors     ;[120,240,50]
;key = ['','h603, no H2','h603, H2','h516, no H2','h516, H2']

IF KEYWORD_SET(moster) THEN $
multiHaloGuo,filesh516[1:2],xrange  = xrange,yrange = yrange,kcorrect = kcorrect,haloids = haloidsh516[*,1:2],psym = psymh516[1:2], color = colorh516[1:2],halomasses = halomasses3,obscolor = obscolor,ctables = ctables[1:2],/things,obssym = obssym,moster = 'moster.stars.z0' ELSE $
multiHaloGuo,filesh516[1:2],xrange  = xrange,yrange = yrange,kcorrect = kcorrect,haloids = haloidsh516[*,1:2],psym = psymh516[1:2], color = colorh516[1:2],halomasses = halomasses3,obscolor = obscolor,ctables = ctables[1:2],/things,obssym = obssym

multiHaloGuo,filesh603[0],  /overplot,                       kcorrect = kcorrect,haloids = haloidsh603[*,0],  psym = psymh603[0],   color = colorh603[0]                           ,obscolor = obscolor,ctables = ctables[0],/bw;,moster = 'moster.stars.z0';,outplot = outplot
multiHaloGuo,filesh603[1:2],/overplot,                       kcorrect = kcorrect,haloids = haloidsh603[*,1:2],psym = psymh603[1:2], color = colorh603[1:2],halomasses = halomasses1,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh986[0],  /overplot,                       kcorrect = kcorrect,haloids = haloidsh986[*,0],  psym = psymh986[0],   color = colorh986[0]                           ,obscolor = obscolor,ctables = ctables[0],/bw;,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh986[1:2],/overplot,                       kcorrect = kcorrect,haloids = haloidsh986[*,1:2],psym = psymh986[1:2], color = colorh986[1:2],halomasses = halomasses2,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh516[0],  /overplot,                       kcorrect = kcorrect,haloids = haloidsh516[*,0],  psym = psymh516[0],   color = colorh516[0]                           ,obscolor = obscolor,ctables = ctables[0],/bw;,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh516[1:2],/overplot,                       kcorrect = kcorrect,haloids = haloidsh516[*,1:2],psym = psymh516[1:2], color = colorh516[1:2],halomasses = halomasses3,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh799[0],  /overplot,                       kcorrect = kcorrect,haloids = haloidsh799[*,0],  psym = psymh799[0],   color = colorh799[0]                           ,obscolor = obscolor,ctables = ctables[0],/bw;,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh799[1:2],/overplot,                       kcorrect = kcorrect,haloids = haloidsh799[*,1:2],psym = psymh799[1:2], color = colorh799[1:2],halomasses = halomasses4,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot

IF KEYWORD_SET(color) THEN loadct,39 ELSE loadct,0
;key = ['no Metals','Metals, no H2', 'H2']
key   = ['no Metals, no ' + textoidl('H_2'),'Metals, no ' + textoidl('H_2'),  textoidl('Metals, H_2')]
key   = ['Primordial','Metals','Metals + ' + textoidl('H_2')]
;key = ['h603, no H2','h603, H2','h516, no H2','h516, H2']
;IF keyword_set(moster) THEN legend,['Moster','THINGS Obs',key],psym = [-3,obssym,psymh603],color = [fgcolor,obscolor,colorh603],/bottom,/right $
;                       ELSE legend,['THINGS Obs',key],psym = [obssym,psymh603],color = [obscolor,colorh603],/bottom,/right 
;IF keyword_set(moster) THEN legend,['Moster','THINGS Obs',key],psym = [-3,obssym,psymh603[2]],color = [fgcolor,obscolor,colorh603[2]],/bottom,/right $
;                       ELSE legend,['THINGS Obs',key[2]],psym = [obssym,psymh603[2]],color = [obscolor,colorh603[2]],/bottom,/right
IF keyword_set(moster) THEN BEGIN 
;    legend,['Moster','THINGS Obs'],psym = [-3,obssym],color = [fgcolor,obscolor],/bottom,/right,box = 0 
    n = 5
ENDIF ELSE BEGIN
;    legend,['THINGS Obs'],psym = [obssym],color = [obscolor],/bottom,/right,box = 0
    n = 4
ENDELSE
lb_keys = strarr(n)
;lb_thicks = intarr(n)
lb_psym = intarr(n) + 3
lb_color = intarr(n) + bgcolor
lb_symsizes = intarr(n) + 1
lb_linestyle = intarr(n)
IF keyword_set(moster) THEN BEGIN 
    l_keys = lb_keys
    l_keys[0:1] = ['Moster','THINGS Obs.']
    l_color = lb_color
    l_color[0:1] = [fgcolor,obscolor]
    l_psym = lb_psym
    l_psym[0:1] = [-3,obssym]
    l_linestyle = lb_linestyle
    legend,l_keys,color = l_color,linestyle = l_linestyle,psym = l_psym,/right,/bottom,box = 0
ENDIF ELSE BEGIN
    l_keys = lb_keys
    l_keys[0] = 'THINGS Obs.'
    l_color = lb_color
    l_color[0] = obscolor
    l_psym = lb_psym
    l_psym[0] = obssym
    l_linestyle = lb_linestyle
    legend,l_keys,color = l_color,linestyle = l_linestyle,psym = l_psym,/right,/bottom,box = 0
ENDELSE
FOR iFile = n - 3, n - 1 DO BEGIN
    loadct,ctables[iFile - (n - 3)]
    l_keys = lb_keys
    l_keys[iFile] = key[iFile - (n - 3)]
;        l_thicks = lb_thicks
;        l_thicks[iFile + 1] = thicks[iFile]
    l_color = lb_color
    l_color[iFile] = colors[iFile - (n - 3)]
    l_psym = lb_psym
    l_psym[iFile] = psym[iFile - (n - 3)]
    l_symsizes = lb_symsizes
    l_symsizes[0]  = 1 
    l_linestyle = lb_linestyle
    legend,l_keys,color = l_color,linestyle = l_linestyle,psym = l_psym,/right,/bottom,box = 0
ENDFOR

;----------------------- Redshift 4 -------------------------------
IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
    stop
    device, filename=outplot + '_guoz4.eps',/COLOR,bits_per_pixel= 8,/times,ysize=12,xsize=18;,/inch
ENDIF ELSE BEGIN
    stop
    window,2
ENDELSE

filesh603 = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.00072.dir/h603.cosmo50cmb.2304g2bwK.00072' ,$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072gs1MbwK/steps/h603.cosmo50cmb.3072gs1MbwK.00072.dir/h603.cosmo50cmb.3072gs1MbwK.00072',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00072.dir/h603.cosmo50cmb.3072g14HBWK.00072']
haloidsh603= [[1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1],$
              [1,2,3,4,5,6,7,8, 9,10,11,12,13,14,15,16,17],$
              [1,2,3,5,6,4,7,8,10,13,11,14,12, 9,15,17,16]]
filesh986 = ['/astro/net/astro-quinn/fabio/REPOSITORY/e11Gals/h986.cosmo50cmb.3072gs1MbwK/h986.cosmo50cmb.3072gs1MbwK.00072/h986.cosmo50cmb.3072gs1MbwK.00072',$
               '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00072.dir/h986.cosmo50cmb.3072g14HBWK.00072']
haloidsh986 = [[1,2,3,4,5,6,7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21],$
               [1,2,3,4,5,7,6,15,11,10,20,13,17,16,19, 9,18,14,22,30,22]]
filesh516 = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00072.dir/h516.cosmo25cmb.3072g1MBWK.00072',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00080.dir/h516.cosmo25cmb.3072g14HBWK.00080']
haloidsh516 = [[1,2,3,4,6],$
               [2,1,4,3,6]]
filesh799 = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g1MBWK/steps/h799.cosmo25cmb.3072g1MBWK.00072.dir/h799.cosmo25cmb.3072g1MBWK.00072',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h799.cosmo25cmb.3072g/h799.cosmo25cmb.3072g14HBWK/steps/h799.cosmo25cmb.3072g14HBWK.00072.dir/h799.cosmo25cmb.3072g14HBWK.00072']
haloidsh799 = [[1,2,3,5,7],$
               [1,2,3,5,8]]

xrange = [8.5,10.5]
yrange = [5,11]
psymh603  = psym     ;[4,4,4]
colorh603 = colors    ;[120,240,50]
psymh986  = psym[1:2]     ;[4,4,4]
colorh986 = colors[1:2]    ;[120,240,50]
psymh516  = psym[1:2] ;[5,5]
colorh516 = colors[1:2];[240,50]
psymh799  = psym[1:2] ;[5,5]
colorh799 = colors[1:2];[240,50]
;key = ['','h603, no H2','h603, H2','h516, no H2','h516, H2']
IF KEYWORD_SET(moster) THEN $
multiHaloGuo,filesh516,    xrange  = xrange,yrange = yrange,kcorrect = kcorrect,haloids = haloidsh516,       psym = psymh516,      color = colorh516,     halomasses = halomasses7,obscolor = obscolor,ctables = ctables[1:2],moster = 'moster.stars.z3.5' ELSE $
multiHaloGuo,filesh516,    xrange  = xrange,yrange = yrange,kcorrect = kcorrect,haloids = haloidsh516,       psym = psymh516,      color = colorh516,     halomasses = halomasses7,obscolor = obscolor,ctables = ctables[1:2]

multiHaloGuo,filesh603[1:2],/overplot,                       kcorrect = kcorrect,haloids = haloidsh603[*,1:2],psym = psymh603[1:2], color = colorh603[1:2],halomasses = halomasses5,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh986,     /overplot,                       kcorrect = kcorrect,haloids = haloidsh986,       psym = psymh986,      color = colorh986,     halomasses = halomasses6,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh516,     /overplot,                       kcorrect = kcorrect,haloids = haloidsh516,       psym = psymh516,      color = colorh516,     halomasses = halomasses7,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh799,     /overplot,                       kcorrect = kcorrect,haloids = haloidsh799,       psym = psymh799,      color = colorh799,     halomasses = halomasses8,obscolor = obscolor,ctables = ctables[1:2];,xrange = xrange,yrange = yrange;,outplot = outplot

;IF KEYWORD_SET(color) THEN loadct,39 ELSE loadct,0
;key = ['no Metals','Metals, no H2', 'H2']
;key = ['h603, no H2','h603, H2','h516, no H2','h516, H2']
;legend,['Moster (z = 3.5)',key[0:2]],psym = [-3,psymh603],color = [fgcolor,colorh603],/bottom,/right
;legend,[key[0:2]],psym = [psymh603],color = [colorh603],/bottom,/right

;----------------------------- Ratio

IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
    device, filename=outplot + '_guoratio.eps',/COLOR,bits_per_pixel= 8,/times,ysize=12,xsize=18;,/inch
ENDIF ELSE BEGIN
    stop
    window,2
ENDELSE

plot, alog10(halomasses1[*,0]),alog10(halomasses1[*,1]),xtitle = textoidl('log(M_{halo}[M') + sunsymbol() + '])',ytitle = textoidl('log(M_{*, Metals + H2}/M_{*, Metals})'),psym = 8,xrange = [9,11.5],yrange = [-0.4,0.4];[-0.4,0.4];,yrange = [-0.6,0.6]
oplot,alog10(halomasses2[*,0]),alog10(halomasses2[*,1]),psym = 8
oplot,alog10(halomasses3[*,0]),alog10(halomasses3[*,1]),psym = 8
oplot,alog10(halomasses4[*,0]),alog10(halomasses4[*,1]),psym = 8
;oplot,alog10(halomasses5[*,0]),alog10(halomasses5[*,1]),psym = 4
;oplot,alog10(halomasses6[*,0]),alog10(halomasses6[*,1]),psym = 4
;oplot,alog10(halomasses7[*,0]),alog10(halomasses7[*,1]),psym = 4
;oplot,alog10(halomasses8[*,0]),alog10(halomasses8[*,1]),psym = 4
oplot,[8,12],[0,0],linestyle = 2
;legend,['h604','h516'],psym = [psym1[0],psym2[0]]
;legend,['z = 0','z = 4'],psym = [8,4],/top,/right

IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
END


PRO twinsanalysis,outplot = outplot2,color = color, kcorrect = kcorrect
;outplot = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/Twins.3072g.guo.324.kcorrect.eps']

;at least 40,000 DM particles
formatplot,outplot = outplot2
IF KEYWORD_SET(outplot2) THEN BEGIN
    fgcolor = 0.1
    if keyword_set(kcorrect) THEN outplot = outplot2 + 'obs' else outplot = outplot2
    device, filename=outplot + '_guoz0.eps',/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
ENDIF ELSE BEGIN
    fgcolor = 255
    window,0
ENDELSE
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL

filesh603 = ['/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.2304g2bwK/steps/h603.cosmo50cmb.2304g2bwK.00328.dir/h603.cosmo50cmb.2304g2bwK.00328',$
             '/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.00324/h603.cosmo50cmb.3072gs1MbwK.00324',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324']
;haloidsh603 = [[1, 2, 3, 6, 17, 12, 15, 16, 25, 44,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1], $
;               [1, 2, 3, 4,  5,  6,  7,  8,  9, 15, 16, 18, 20, 21, 23, 26, 30, 31, 33, 34, 35, 36, 38, 39, 40, 41, 42, 43, 44, 52],$
;               [1, 2, 3, 4,  5,  7,  6,  8,  9, 14, 16, 18, 19, 22, 24, 46, 31, 30, 33, 36, 34, 35, 39, 40, 41, 38, 43, 44, 49, 53]]
haloidsh603 = [[1, 2, 3, 6, 17, 12, 15, 16, 25, 44,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1], $
               [1, 2, 3, 4,  5,  6,  7,  8,  9, 15, 16, 18, 20, 21, 23, 26, 30, 31, 33, 34],$
               [1, 2, 3, 4,  5,  7,  6,  8,  9, 14, 16, 18, 19, 22, 24, 46, 31, 30, 33, 36]]
;52: 20315 -- valid to 425
;76: 13845, 53:20714
;haloidsh603 = [['1', '2', '3', '6', '17', '12', '15', '16', '25', '44',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1',  '1'], $
;               ['1', '2', '3', '4',  '5',  '6',  '7',  '8',  '9', '15', '16', '18', '20', '21', '23', '26', '30', '31', '33', '34', '35', '36', '38', '39', '40', '41', '42', '43', '44', '51', '52'],$
;               ['1', '2', '3', '4',  '5',  '7',  '6',  '8',  '9', '14', '16', '18', '19', '22', '24', '46', '31', '30', '33', '36', '34', '35', '39', '40', '41', '38', '43', '44', '49', '76', '53']]
filesh516 = ['/astro/net/nbody1/abrooks/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.00324/h516.cosmo25cmb.3072g1MBWK.00324',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00320.dir/h516.cosmo25cmb.3072g14HBWK.00320']
;x = matchHalos(filesh516 + '.amiga.stat')
;haloidsh516 = [[1,2,3,4,6,9,10,15,18,28,29],$
;               [1,2,3,4,6,9,13,15,18,28,23]]
haloidsh516 = [[1,2,3,4,6,9,10],$
               [1,2,3,4,6,9,13]]
;haloids = [[1,2,3,4,5,6,15,16,18,20],$
;           [1,2,3,4,5,7,14,16,18,19]]

xrange = [8.5,12]
yrange = [5,11]

n = N_ELEMENTS(filesh603)
IF KEYWORD_SET(color) THEN BEGIN
    loadct,39
    obscolor = 0 ;fgcolor
    if colors[0] eq 1 then  colors  = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(psym) THEN psym = fltarr(n) + 4
ENDIF ELSE BEGIN
    loadct,0    
    obscolor = 100
    colors = (fltarr(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
    IF NOT KEYWORD_SET(psym) THEN  psym = (findgen(n)+2)*2
ENDELSE
psymh603  = psym      ;[4,4,4]
colorh603 = colors     ;[120,240,50]
psymh516  = psym[1:2] ;[5,5]
colorh516 = colors[1:2];[240,50]
;key = ['','h603, no H2','h603, H2','h516, no H2','h516, H2']
multiHaloGuo,filesh603[0],  xrange  = xrange,yrange = yrange,kcorrect = kcorrect,haloids = haloidsh603[*,0],  psym = psymh603[0],   color = colorh603[0]                           ,obscolor = obscolor,/bw,/things,moster = 'moster.stars.z0.5';,outplot = outplot
multiHaloGuo,filesh603[1:2],/overplot,                       kcorrect = kcorrect,haloids = haloidsh603[*,1:2],psym = psymh603[1:2], color = colorh603[1:2],halomasses = halomasses1,obscolor = obscolor,/bw;,xrange = xrange,yrange = yrange;,outplot = outplot
multiHaloGuo,filesh516,     /overplot,                       kcorrect = kcorrect,haloids = haloidsh516,       psym = psymh516,      color = colorh516,     halomasses = halomasses2,obscolor = obscolor,/bw;,xrange = xrange,yrange = yrange;,outplot = outplot

IF KEYWORD_SET(color) THEN loadct,39 ELSE loadct,0
key = ['no Metals','Metals, no H2', 'H2']
;key = ['h603, no H2','h603, H2','h516, no H2','h516, H2']
legend,['Moster','THINGS Obs.',key],psym = [-3,2,psymh603],color = [fgcolor,obscolor,colorh603],/bottom,/right


IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
    device, filename=outplot + '_guoz4.eps',/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
ENDIF ELSE BEGIN
    stop
    window,2
ENDELSE

;----------------------- Redshift 4 -------------------------------
filesh603 = ['/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.00072/h603.cosmo50cmb.3072gs1MbwK.00072',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00072.dir/h603.cosmo50cmb.3072g14HBWK.00072']
haloidsh603 = [['1','2','3','6','4','5','7','8','14','9','11','13','10','12','15','17','16','23','19','18','25','28','20','26','21','24','27','26','22','29','30','40','33','37','38','32'], $
               ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','22','23','24','25','26','27','28','29','30','31','32','33','35','36','37','39']]
;40:42549
;41:44518, 42:35276
;,'21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','54','57'], $
;,'21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','39','40','41','54','57']]
filesh516 = ['/astro/net/nbody1/abrooks/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.00072/h516.cosmo25cmb.3072g1MBWK.00072',$
             '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00080.dir/h516.cosmo25cmb.3072g14HBWK.00080']
;haloidsh516 = [[1,2,3,4,6,11],$
;               [2,1,4,3,6, 7]]
haloidsh516 = [[1,2,3,4,6],$
               [2,1,4,3,6]]
xrange = [8.5,10.5]
yrange = [5,11]
n = N_ELEMENTS(filesh603)
;IF KEYWORD_SET(color) THEN BEGIN
;    loadct,39
;    obscolor = 0 ;fgcolor
;    if colors[0] eq 1 then  colors  = (findgen(n) + 1)*240/n else colors = colors
;    IF NOT KEYWORD_SET(psym) THEN psym = fltarr(n) + 4
;ENDIF ELSE BEGIN
;    loadct,0    
;    obscolor = 100
;    colors = (fltarr(n) + 1)*fgcolor ;(findgen(n) + 1)*10.0 + 5.0;  fltarr(N_ELEMENTS(broadband)) + 5
;    IF NOT KEYWORD_SET(psym) THEN  psym = (findgen(n)+2)*2
;ENDELSE
psymh603  = psym[1:2]     ;[4,4,4]
colorh603 = colors[1:2]    ;[120,240,50]
psymh516  = psym[1:2] ;[5,5]
colorh516 = colors[1:2];[240,50]
;key = ['','h603, no H2','h603, H2','h516, no H2','h516, H2']
multiHaloGuo,filesh603,  xrange  = xrange,yrange = yrange,kcorrect = kcorrect,haloids = haloidsh603,psym = psymh603,color = colorh603,halomasses = halomasses3,obscolor = obscolor,/bw;,outplot = outplot
multiHaloGuo,filesh516,     /overplot,                    kcorrect = kcorrect,haloids = haloidsh516,psym = psymh516,color = colorh516,halomasses = halomasses4,obscolor = obscolor,/bw;,xrange = xrange,yrange = yrange;,outplot = outplot

IF KEYWORD_SET(color) THEN loadct,39 ELSE loadct,0
;key = ['no Metals','Metals, no H2', 'H2']
;key = ['h603, no H2','h603, H2','h516, no H2','h516, H2']
;legend,['Moster (z = 0.5)','THINGS Obs (z = 0)',key[1:2]],psym = [-3,2,psymh603],color = [fgcolor,obscolor,colorh603],/bottom,/right
legend,[key[1:2]],psym = [psymh603],color = [colorh603],/bottom,/right
IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
    device, filename=outplot + '_guoratio.eps',/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
ENDIF ELSE BEGIN
    stop
    window,1
ENDELSE

plot, alog10(halomasses1[*,0]),alog10(halomasses1[*,1]),xtitle = 'log(Mass Halo)',ytitle = textoidl('log(M_{*, H2}/M_{* Metals})'),psym = 8,xrange = [8.5,12],yrange = [-1,1];[-0.4,0.4];,yrange = [-0.6,0.6]
oplot,alog10(halomasses2[*,0]),alog10(halomasses2[*,1]),psym = 8
oplot,alog10(halomasses3[*,0]),alog10(halomasses3[*,1]),psym = 4
oplot,alog10(halomasses4[*,0]),alog10(halomasses4[*,1]),psym = 4
oplot,[8,12],[0,0],linestyle = 2
;legend,['h604','h516'],psym = [psym1[0],psym2[0]]
;legend,['z = 0','z = 4'],psym = [8,4]

IF KEYWORD_SET(outplot) THEN device,/close ELSE stop
END

;twinsanalysis_ZSTEPS,outplot = '~/plots/h603/stellarmass_halomass'
PRO twinsanalysis_zsteps,outplot = outplot2
formatplot,outplot = outplot2
IF KEYWORD_SET(outplot2) THEN BEGIN
    fgcolor = 0.1
    if keyword_set(kcorrect) THEN outplot = outplot2 + 'obs' else outplot = outplot2
    device, filename=outplot + '_guozall.eps',/COLOR,bits_per_pixel= 8,/times,ysize=5,xsize=7,/inch
ENDIF ELSE BEGIN
    fgcolor = 255
    window,0
ENDELSE

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/'
base = 'h603.cosmo50cmb.3072g14HBWK'
steps =   ['00072','00084','00096','00108','00120','00132','00144','00156','00168','00180','00192','00204','00216','00228','00240','00252','00264','00276','00288','00300','00312','00324']
z =       ['3.436','2.999','2.654','2.372','2.138','1.939','1.766','1.616','1.483','1.364','1.257','1.160','1.072','0.991','0.916','0.846','0.782','0.722','0.665','0.612','0.562','0.515']
tags =    ['3.5'  ,    ' ',    ' ',    ' ',    ' ',    '2',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    '1',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',  '0.5']
haloids = [[    2],[    2],[    2],[    2],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1]]
offset = [-0.1,0.1]
imf = 1
nsteps = N_ELEMENTS(steps);last/dt
files = strarr(nsteps)
key = strarr(nsteps)
color = fltarr(nsteps) + fgcolor
psym = fltarr(nsteps) - 3
psym[where(tags ne ' ')] = -4
symsize = fltarr(nsteps) + 1
xrange = [8.5,12]
yrange = [5,11]

FOR i = 0, nsteps - 1 DO BEGIN 
    step = steps[i]
    files[i] = prefix + base + '/steps/' + base + '.' + step + '.dir/' + base + '.' + step
    key[i] = step
ENDFOR
multiHaloGuo,files,yrange = yrange,xrange = xrange,psym = psym,imf = imf,haloids = haloids,tags = tags,offset = offset,mainhalos = mainhalosh603,/bw,color = color;,/kcorrect;,key = key,
oplot,alog10(mainhalosh603[*,0]),alog10(mainhalosh603[*,1]);,psym = psym[0]

prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/'
base = 'h516.cosmo25cmb.3072g14HBWK'
steps =   ['00024','00032','00080','00088','00092','00104','00128','00132','00144','00148','00156','00168','00180','00192','00204','00216','00228','00240','00252','00260','00272','00284','00308','00320','00332','00344','00356','00368','00380','00392','00404','00420','00432','00436','00444','00456','00468','00480','00492','00504','00512']
tags =    [' '    ,    ' ',   '3.5',    ' ',    ' ',    ' ',    '2',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    '1',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',   '0.5',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    ' ',    '0']
haloids = [[    1],[    1],[    1],[     1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[     1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1],[    1]]
offset = [-0.1,0.2]
imf = 1
nsteps = N_ELEMENTS(steps);last/dt
files = strarr(nsteps)
key = strarr(nsteps)
color = fltarr(nsteps)
psym = fltarr(nsteps) - 3
psym[where(tags ne ' ')] = -6
symsize = fltarr(nsteps) + 1
FOR i = 0, nsteps - 1 DO BEGIN 
    step = steps[i]
    files[i] = prefix + base + '/steps/' + base + '.' + step + '.dir/' + base + '.' + step
    key[i] = step
ENDFOR
multiHaloGuo,files,yrange = yrange,xrange = xrange,psym = psym,imf = imf,haloids = haloids,tags = tags,/overplot,mainhalos = mainhalosh516,/bw,color = color;,/kcorrect;,key = key,
oplot,alog10(mainhalosh516[*,0]),alog10(mainhalosh516[*,1]);,psym = psym[0]
IF KEYWORD_SET(outplot) THEN BEGIN
    device,/close
ENDIF ELSE BEGIN
    stop
ENDELSE
END
