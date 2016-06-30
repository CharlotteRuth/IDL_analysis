
pro sfh,files,pfiles,keys = keys,colors = colors,linestyle = linestyle,thick = thick,ctables = ctables,yrange = yrange,xrange = xrange,binsize = binsize,cumlative = cumlative,label = label, redshift = redshift,outplot = outplot,formatthick = formatthick,bulge = bulge,_EXTRA=_extra
n = N_ELEMENTS(files)

!Y.STYLE = 1
!X.STYLE = 1
formatplot,outplot = outplot,thick = formatthick
IF KEYWORD_SET(outplot) THEN BEGIN
    fgcolor = 0
    bgcolor = 255
    device,filename=outplot + '_sfh.eps',/color,bits_per_pixel= 8,/times,xsize = 18,ysize= 12
ENDIF ELSE BEGIN
    fgcolor = 255
    bgcolor = 0
    window,0
ENDELSE

IF KEYWORD_SET(colors) THEN BEGIN
    loadct,39
    if NOT keyword_set(ctables) then ctables = intarr(n) + 39
    if colors[0] eq 1 then  colors = (findgen(n) + 1)*240/n else colors = colors
    IF NOT KEYWORD_SET(thick) THEN thick = fltarr(n) + 2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = fltarr(n) 
ENDIF ELSE BEGIN
    loadct,0    
    colors = fltarr(n) + fgcolor
    if NOT keyword_set(ctables) then ctables = fltarr(n)
    IF NOT KEYWORD_SET(thick) THEN thick = findgen(n)*2
    IF NOT KEYWORD_SET(linestyle) THEN linestyle = REVERSE(findgen(n)*2   )
ENDELSE
IF NOT KEYWORD_SET(binsize) THEN binsize = 5e7
IF KEYWORD_SET(xrange) THEN maxt =   xrange[1]*1e9
print,linestyle

lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)

FOR i = 0, n - 1 DO BEGIN
    loadct,ctables[i]
    units = tipsyunits(pfiles[i])
;    stop
   rtipsy,files[i],h,g,d,s
   IF i eq 0 THEN  sfr,s,massunit = units.massunit,timeunit = units.timeunit,binsize=binsize,thick = thick[i],linestyle = linestyle[i],cumlative = cumlative,sarray = sarray,tarray = tarray,maxt = maxt,title = label,yrange =yrange,_EXTRA=_extra,/nodata,redshift = redshift
   sfr,s,                massunit = units.massunit,timeunit = units.timeunit,binsize=binsize,thick = thick[i],linestyle = linestyle[i],cumlative = cumlative,sarray = sarray,tarray = tarray,maxt = maxt,color = colors[i],/overplot
   ind = where(tarray/1e9 gt 10)
   print,'Mean: ',mean(sarray(ind))
   print,'STDEV: ',stdev(sarray(ind))

   IF keyword_set(bulge) THEN BEGIN
       readarr,files[i] + '.decomp',h,deco,part = 'star',type = 'long',/ascii
       bulge = s[where(deco EQ 3 OR deco EQ 5)]
       sfr,bulge,massunit = units.massunit,timeunit = units.timeunit,binsize=binsize,thick = thick[i] - 4,linestyle = linestyle[i],cumlative = cumlative,sarray = sarray,tarray = tarray,maxt = maxt,color = colors[i],/overplot       
       print,n_elements(where(bulge.tform*units.timeunit/1e9 GE 6 AND bulge.tform*units.timeunit/1e9 LE 8)),n_elements(bulge),n_elements(where(bulge.tform*units.timeunit/1e9 GE 6 AND bulge.tform*units.timeunit/1e9 LE 8))/float(n_elements(bulge))
       print,n_elements(where(bulge.tform*units.timeunit/1e9 GE 6 AND bulge.tform*units.timeunit/1e9 LE 8)),n_elements(where(s.tform*units.timeunit/1e9 GE 6 AND s.tform*units.timeunit/1e9 LE 8)),n_elements(where(bulge.tform*units.timeunit/1e9 GE 6 AND bulge.tform*units.timeunit/1e9 LE 8))/float(n_elements(where(s.tform*units.timeunit/1e9 GE 6 AND s.tform*units.timeunit/1e9 LE 8)))
   ENDIF

 ;  IF KEYWORD_SET(keys) THEN BEGIN
 ;      l_keys = lb_keys
 ;      l_keys[i] = keys[i]
 ;      l_thicks = lb_thicks
 ;      l_thicks[i] = thick[i]
 ;      l_color = lb_color
 ;      l_color[i] = colors[i]
 ;      l_linestyle = lb_linestyle
 ;      l_linestyle[i] = linestyle[i]
 ;      legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/right,/top,box = 0
 ;  ENDIF
ENDFOR
IF KEYWORD_SET(keys) THEN legend,keys,color = colors,linestyle = linestyle,thick = thick,ctables = ctables,/right,/top,box = 0
if keyword_set(outplot) then device,/close else stop
end

pro sfh_master

;dir0 = '/astro/net/nbody1/jillian/h516/'
;file0 = 'h516.cosmo25cmb.1536g1MBWK.00512'
dir0 = '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g3HBWK/steps/old_fb/h603.cosmo50cmb.1536g3HBWK.00144.dir/'
file0 = 'h603.cosmo50cmb.1536g3BWK.00144'

;dir1 = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00512.dir/'
;file1 = 'h516.cosmo25cmb.1536g2HBWK.00512'
dir1 = '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g2bwK.BUG/h603.cosmo50cmb.2304g2bwK.00144/'
file1 = 'h603.cosmo50cmb.2304g2bwK.00144'

;dir2 = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00513.dir/'
;file2 = 'h516.cosmo25cmb.1536g3HBWK.00513'
dir2 = '/astro/net/scratch2/fabio/REPOSITORY/e11Gals/h603.cosmo50cmb.2304g5bwK.BUG/h603.cosmo50cmb.2304g5bwK.00144.dir/'
file2 = 'h603.cosmo50cmb.2304g5bwK.00144'

dir3 = '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g3HBWK/steps/h603.cosmo50cmb.1536g3HBWK.00180.dir/'
file3 = 'h603.cosmo50cmb.1536g3HBWK.00180'

;pfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g3HBWK/h516.cosmo25cmb.1536g3HBWK.param'
pfile = '/astro/net/scratch2/christensen/MolecH/Cosmo/h603.cosmo50cmb.1536g3HBWK/h603.cosmo50cmb.1536g3HBWK.param'

;key = ['No H2','H2','H2 SF']
key = ['H2','low threshhold','highthreshold','new H2']
colors = [140,240,50,100]
linestyle = [0,0,0,0]
yrange = [0,5]
;sfh,files,pfiles,key = key,color = color,linestyle = linestyle,outplot = outplot,color = color,thicks = thicks,_EXTRA=_extra

;-------------- h516------------
prefix = '/Users/christensen/Data/MolecH/Cosmo/'
IF KEYWORD_SET(vetinari) THEN prefix = '/Users/christensen/Data/MolecH/Cosmo/' ELSE prefix = '/astro/users/christensen/Scratch/MolecH/Cosmo/' 
files = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir/h516.cosmo25cmb.3072g1MBWK.00492.halo.1',$
         'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1'] 
key = ['DnoH2','DH2']
IF KEYWORD_SET(vetinari) THEN outplot = '~/Figures/h516.cosmo25cmb.paper.sfh0.3.eps' ELSE outplot = '~/plots/h516.cosmo25cmb.paper.sfh.eps'
pfiles = prefix + ['h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param', $
                   'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']


;--------- Twins -------------------
files = ['/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.00324/h603.cosmo50cmb.3072gs1MbwK.00324.1.std',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/steps/h603.cosmo50cmb.3072g14HBWK.00324.dir/h603.cosmo50cmb.3072g14HBWK.00324.halo.1']
key = ['h603, no H2','h603, H2']
pfiles =  ['/astro/net/nbody1/abrooks/h603.cosmo50cmb.3072gs1MbwK/h603.cosmo50cmb.3072gs1MbwK.param',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g14HBWK/h603.cosmo50cmb.3072g14HBWK.param']

outplot = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h603.cosmo50cmb.3072g/h603.cosmo50cmb.3072g.CumSFH.eps'
sfh,files,pfiles,key = key,outplot = outplot,/color,/cumlative

files = ['/astro/net/nbody1/abrooks/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.00512/h516.cosmo25cmb.3072g1MBWK.00512.1.std',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/h516.cosmo25cmb.3072g14HBWK.00512.halo.1']
key = ['h516, no H2','h516, H2']
pfiles =  ['/astro/net/nbody1/abrooks/h516.cosmo25cmb.3072g1MBWK/h516.cosmo25cmb.3072g1MBWK.param',$
         '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param']
outplot = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g.CumSFH.eps'

sfh,files,pfiles,key = key,outplot = outplot,/color,/cumlative;,yrange = [0,0.3]

end
