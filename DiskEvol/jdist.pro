pro jdist, prefix=prefix, filename, h, hs, hg, d, s, g, fdisk, dmass, smass, gmass, halo = halo, rmax=rmax, all=all, obs=obs, disk=disk

boxsize = strmid(filename,5,7)
if boxsize eq 'cosmo25' then lunit=25000.
if boxsize eq 'cosmo50' then lunit=50000.

if keyword_set(prefix) then filename=prefix+filename
IF NOT keyword_set(halo) THEN halo = 1
halo_str = strtrim(halo,2)

IF NOT (FILE_TEST(filename + ".halo."+halo_str+".std")) THEN BEGIN
    rtipsy, filename, h,g,d,s
    grp = read_lon_array(filename+'.amiga.grp')
    if exist(filename+'.cmp') then cmp = read_lon_array(filename+'.cmp')
    rtipsy, filename+'.amiga.gtp', h1,g1,d1,s1
    readarr, filename + '.HI',h,hi,/ascii,part = 'gas'
    IF (FILE_TEST(filename + ".H2")) THEN  readarr,filename+".H2",h,H2,/ascii,part = 'gas' ELSE H2 = HI * 0
;Get main halo
    ggrp = where(grp[0:h.ngas-1] eq halo, ng)
    dgrp = where(grp[h.ngas:h.ndark+h.ngas-1] eq halo, nd)
    sgrp = where(grp[h.ngas+h.ndark:h.n-1] eq halo, ns)
    g = g[ggrp]
    HI = HI[ggrp]
    H2 = 2.0*H2[ggrp]
    d = d[dgrp]
;IF KEYWORD_SET(all) eq 0 THEN BEGIN
;ntot = ng+nd+ns
    
    IF KEYWORD_SET(disk) THEN begin
;if exist(filename+'.cmp') then begin
        cmp = cmp[h.ngas+h.ndark:h.n-1]
;For bulge + disk
;IF KEYWORD_SET(obs) THEN scmp = where(cmp eq 1 or cmp eq 3)
;For disk alone
        scmp = where(cmp eq 1)
;bulgemass = total(s[where(cmp eq 3)].mass)
        IF h.nstars NE 0 THEN s = s[scmp]
;ENDIF
;endif
    ENDIF ELSE s = s[sgrp]
;ns = n_elements(s)
    
;Center on main halo
    g.x = g.x-s1[0].x
    g.y = g.y-s1[0].y
    g.z = g.z-s1[0].z
    d.x = d.x-s1[0].x
    d.y = d.y-s1[0].y
    d.z = d.z-s1[0].z
    s.x = s.x-s1[0].x
    s.y = s.y-s1[0].y
    s.z = s.z-s1[0].z
    g.vx = g.vx-s1[0].vx
    g.vy = g.vy-s1[0].vy
    g.vz = g.vz-s1[0].vz
    d.vx = d.vx-s1[0].vx
    d.vy = d.vy-s1[0].vy
    d.vz = d.vz-s1[0].vz
    s.vx = s.vx-s1[0].vx
    s.vy = s.vy-s1[0].vy
    s.vz = s.vz-s1[0].vz
;Align, just for the hell of it
    align, g,d,s, 1.0/lunit
ENDIF ELSE BEGIN
    rtipsy, filename + ".halo." + halo_str + ".std", head,g,d,s
    if exist(filename+'.halo.' + halo_str + '.cmp') then cmp = read_lon_array(filename+'.cmp')
    readarr, filename + '.halo.' + halo_str + '.HI',head,hi,/ascii,part = 'gas'
    readarr, filename + '.halo.' + halo_str + '.HeI',head,hei,/ascii,part = 'gas'
    IF (FILE_TEST(filename + ".halo." + halo_str + ".H2")) THEN  readarr,filename+".H2",head,H2,/ascii,part = 'gas' ELSE H2 = HI * 0
ENDELSE


IF NOT KEYWORD_SET(all) THEN BEGIN
    IF head.nstar NE 0 THEN rads = sqrt(s.x*s.x+s.y*s.y+s.z*s.z)
    radg = sqrt(g.x*g.x+g.y*g.y+g.z*g.z)
    radd = sqrt(d.x*d.x+d.y*d.y+d.z*d.z)
    if not keyword_set(rmax) then rmax = max(rads)
    IF head.nstar NE 0 THEN inds = where(rads lt rmax)
    indd = where(radd lt rmax)
;    indg = where(g.tempg lt 20000.)
    IF head.nstar NE 0 THEN smass = total(s[inds].mass) ELSE smass = 0
;    gmass = total(g[indg].mass)
    gmass = total(g.mass*[HI + H2*2 + HeI])
    dmass = total(d[indd].mass)
    fdisk = (smass+gmass)/dmass 
    fdisks = smass/dmass
;    fbar = 0.2121    ;(omega_b/(omega_m-omega_b)) = 0.042/(0.24-0.042) ;This was the original formulation but Alyson sent me the below correction
    fbar = 0.175      ;(omega_b/omega_m) = 0.042/0.24
ENDIF ELSE BEGIN
    IF head.nstar NE 0 THEN smass = total(s.mass) ELSE smass = 0
    gmass = total(g.mass)
    dmass = total(d.mass)
    fdisk = (smass+gmass)/dmass
    fdisks = smass/dmass
ENDELSE



angmom, d, jvec, lvec, jmag, lmag
;Following Sharma & Steinmetz: they align (independently) the dm and the baryons.
; I have aligned both according to the gas
;They use the z component of the j vector only
;They use only the positive valued j 
; BUT after alignment, are any of my j_z values negative?  <--check
ind = where(jvec[*,2] ge 0.)
d = d[ind]
;h = histogram(jvec[ind,2]/mean(jvec[ind,2]), binsize=0.05, reverse_indices=ri)
h = histogram(jvec[ind,2]/mean(jvec[ind,2]), binsize=0.05, reverse_indices=ri)
;h = histogram(jmag/mean(jmag), binsize=0.05, reverse_indices=ri)
; get the mass in each bin
dmass = fltarr(n_elements(h))
for j=0L,n_elements(h)-1 do if ri[j+1] gt ri[j] then dmass[j] = total(d[ri[ri[j]:ri[j+1]-1]].mass)
;Gas
angmom, g, jvecg, lvec, jmagg, ltot
ind = where(jvecg[*,2] ge 0.)
g = g[ind]
hg = histogram(jvecg[ind,2]/mean(jvecg[ind,2]), binsize=0.05, reverse_indices=ri)
;hg = histogram(jmagg/mean(jmagg), binsize=0.05, reverse_indices=ri)
gmass = fltarr(n_elements(hg))
;Hmass = g.mass*(HI + H2)*1.4
Hmass = g.mass*(HI + HeI + H2)
;for j=0,n_elements(hg)-1 do if ri[j+1] gt ri[j] then gmass[j] = total(g[ri[ri[j]:ri[j+1]-1]].mass)
for j=0,n_elements(hg)-1 do if ri[j+1] gt ri[j] then gmass[j] = total(Hmass[ri[ri[j]:ri[j+1]-1]])

;Stars
IF head.nstar NE 0 THEN BEGIN
    angmom, s, jvecs, lvec, jmags, ltot
    ind = where(jvecs[*,2] ge 0.)
    s = s[ind]
    hs = histogram(jvecs[ind,2]/mean(jvecs[ind,2]), binsize=0.05, reverse_indices=ri)
;hs = histogram(jvecs[ind,2]/mean(jmags), binsize=0.05, reverse_indices=ri)
;hs = histogram(jmags/mean(jmags), binsize=0.05, reverse_indices=ri)
    smass = fltarr(n_elements(hs))
    for j=0,n_elements(hs)-1 do if ri[j+1] gt ri[j] then smass[j] = total(s[ri[ri[j]:ri[j+1]-1]].mass)
;sjjtot = fltarr(n_elements(hs))
;for j=0,n_elements(hs)-1 do if ri[j+1] gt ri[j] then sjjtot[j] = total(jvecs[[ri[ri[j]:ri[j+1]-1]],2])/total(jvecs[ind,2])
;oplot, indgen(n_elements(hs))*0.05, smass/(total(s.mass)*fdisk), linestyle=1
ENDIF ELSE smass = fltarr(n_elements(hg))

end

PRO jdist_plot,outplot = outplot
spawn,'hostname',hostname
IF hostname EQ 'ozma' THEN prefix = '/home/christensen/Storage1/UW/MolecH/Cosmo/' ELSE prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'

dir = ['/astro/store/nbody3/fabio/h986/768.g2',$
       prefix + 'h986.cosmo50cmb.3072g/h986.cosmo50cmb.3072g14HBWK/steps/h986.cosmo50cmb.3072g14HBWK.00512.dir']
filenames = ['h986.cosmo50cmb.768g2.00512',$
             'h986.cosmo50cmb.3072g14HBWK.00512']
halos = [17,1]
cameras = [15,15]               ;45
center = [[0,0],$
          [0,0]]
ctables = [39,39]

cd,dir[1]
rmax = opticalRadii(filename = filenames[1],center = center[*,1],extno = cameras[1],verbose = verbose)

formatplot,outplot = outplot,/thick
IF KEYWORD_SET(outplot) THEN fgcolor = 0 ELSE fgcolor = 255
IF KEYWORD_SET(outplot) THEN bgcolor = 255 ELSE bgcolor = 0
IF keyword_set(outplot) THEN device, filename=outplot + '_jjtot.eps', /encapsulated, /color,/times,ysize = ysize,xsize = xsize,bits_per_pixel = 8  $
ELSE window,0                   ;,ysize = ysize,xsize = xsize
;s =indgen(n_elements(h))*0.1
;xsi = 1.0-1.25*(1.0-(1.25-1.0)*alog(1.25/(1.25-1.0)))
;ps = xsi*1.25*(1.25-1.0)/(xsi*s+1.25-1.0)^2.
loadct,39

n = n_elements(dir)
lb_keys = strarr(n)
lb_thicks = intarr(n)
lb_color = intarr(n) + bgcolor
lb_linestyle = intarr(n)  

FOR i = 1, 1 DO BEGIN
    cd,dir[i]
    loadct,ctables[i]
 ;   print,filenames[i],rmax
;posx=[.15,.95] 
    ytickname=[' ','0.2',' ','0.6',' ','1.0'] 
;posy=[.25,.95]
    jdist, prefix=dir[i] + '/', filenames[i], h, hs, hg, d, s, g, fdisk, dmass, smass, gmass, halo = halos[i],rmax = rmax
    if i eq 1 then plot, indgen(n_elements(h))*0.05, dmass/max(dmass), xrange=[0,2.5], yrange=[0,1], xtickname=xtickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)',/nodata,title = label ;, charthick=2
    if n_elements(hs) lt n_elements(hg) then b = n_elements(hg) else b = n_elements(hs)
    if max(smass) NE 0 THEN BEGIN
;        oplot, indgen(b)*0.05, (gmass/max(gmass)+smass/max(smass))*(fdisk/.2121), color = 30
        oplot, indgen(b)*0.05, (smass/max(smass))*(fdisk/.2121), color = 30 
        polyfill,[0,indgen(b)*0.05,2.5],[0,smass/max(smass)*(fdisk/.2121),0],color = 30
    ENDIF ELSE BEGIN
        oplot, indgen(b)*0.05, gmass/max(gmass)*(fdisk/.2121), color = 30
        polyfill,[0,indgen(b)*0.05,2.5],[0,gmass/max(gmass)*(fdisk/.2121),0],color = 30
    ENDELSE
    oplot,indgen(n_elements(h))*0.05,dmass/max(dmass), color = 130
    legend,['Stars','Dark Matter'],color = [30,130],linestyle = [0,0],/right
;    IF KEYWORD_SET(keys) THEN BEGIN
    IF 0 THEN BEGIN
        l_keys = lb_keys
        l_keys[i] = keys[i]
        l_thicks = lb_thicks
        l_thicks[i] = thicks[i]
        l_color = lb_color
        l_color[i] = colors[i]
        l_linestyle = lb_linestyle
        l_linestyle[i] = 0      ;linestyle[i]
        legend,l_keys,color = l_color,linestyle = l_linestyle,thick = l_thicks,/right,box = 0
    ENDIF
    stop
ENDFOR
IF keyword_set(outplot) THEN device, /close
END
