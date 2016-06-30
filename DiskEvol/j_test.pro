PRO j_test, sim,files,center = center, extno = extno,rmax = rmax

;To be run at UW, using that filesystem
;!p.font=0
;!p.thick=4
;!x.thick=3
;!y.thick=3
aspect_ratio=1.5 ;rectangle
xsize=4
ysize=xsize*aspect_ratio
loadct, 0;13
IF NOT KEYWORD_SET(center) THEN center = [0,0] 
IF NOT KEYWORD_SET(camera) THEN extno = 14 

center = [-8,72]
extno = 38
;sim = 'e10Gals/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.00512/'
;files = 'h516.cosmo25cmb.3072g14HBWK.00512'

;sim = 'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00512.dir/'
;files = 'h516.cosmo25cmb.3072g14HBWK.00512'

;prefix = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
cd,sim ;prefix + sim 
formatplot,/outplot
device, filename=files + '_jjtot.eps', /encapsulated, /color, /inches,/times,ysize = 5,xsize = 8,bits_per_pixel = 8

;center = [-8,72]
;extno = 16
IF NOT KEYWORD_SET(rmax) THEN rmax =  opticalRadii(filename = files,center = center, extno = extno)
;rmax = 8
print,rmax
 ;measured from checkfit -- 
;the point where the surface brightness of the stars drops below 25 mag/arsec^2
;e.g., the "observable" radius of the disk
; extend this out further for dwarfs (3 & 2 kpc in stars for h516 & h799) b/c of large gas disk

jdist, prefix=sim + '/', files, h, hs, hg, d, s, g, fdisk, dmass, smass, gmass, rmax=rmax, /obs 
posx=[.15,.95] 
ytickname=[' ','0.2',' ','0.6',' ','1.0'] 
posy=[.25,.95]

;FOR a generic DM halo:
s =indgen(n_elements(h))*0.1
xsi = 1.0-1.25*(1.0-(1.25-1.0)*alog(1.25/(1.25-1.0)))
ps = xsi*1.25*(1.25-1.0)/(xsi*s+1.25-1.0)^2.
;plot, s, ps, xrange=[0,4], yrange=[0,1], pos=[posx[0],posy[0],posx[1],posy[1]], xtickname=xtickname, ytickname=ytickname, charthick=2, charsize=2


plot, indgen(n_elements(h))*0.05, dmass/max(dmass), xrange=[0,4], yrange=[0,1], pos=[posx[0],posy[0],posx[1],posy[1]], xtickname=xtickname, ytickname=ytickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)';, charthick=2
if n_elements(hs) lt n_elements(hg) then b = n_elements(hg) else b = n_elements(hs)
oplot, indgen(b)*0.05, (gmass/max(gmass)+smass/max(smass))*(fdisk/.2121), linestyle=2
legend, ['DM', 'Observable'], linestyle=[0,2], /top, /right ;charthick = 2
;xyouts, 1.8, -0.2, textoidl('s = j/j_{tot}'), charthick=3, charsize=1.5
;xyouts, -.5, 0.45, 'P(s)', charthick=3, orientation=90, charsize=1.5
device, /close


device, filename=files + '_jjtot_stars.eps', /encapsulated, /color, /inches,/times,ysize = 5,xsize = 8,bits_per_pixel = 8
plot, indgen(n_elements(h))*0.05, dmass/max(dmass), xrange=[0,4], yrange=[0,1], pos=[posx[0],posy[0],posx[1],posy[1]], xtickname=xtickname, ytickname=ytickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)';, charthick=2                                                                                                                                               
if n_elements(hs) lt n_elements(hg) then b = n_elements(hg) else b = n_elements(hs)
oplot, indgen(b)*0.05,(smass/max(smass))*(fdisk/.2121), linestyle=2
legend, ['DM', 'Stars'], linestyle=[0,2], /top, /right ;charthick = 2                                                                          
;xyouts, 1.8, -0.2, textoidl('s = j/j_{tot}'), charthick=3,
;charsize=1.5                                                                                                                
;xyouts, -.5, 0.45, 'P(s)', charthick=3, orientation=90, charsize=1.5                                                                                                                   
device, /close

device, filename=files + '_jjtot_gas.eps', /encapsulated, /color, /inches,/times,ysize = 5,xsize = 8,bits_per_pixel = 8
plot, indgen(n_elements(h))*0.05, dmass/max(dmass), xrange=[0,4], yrange=[0,1], pos=[posx[0],posy[0],posx[1],posy[1]], xtickname=xtickname, ytickname=ytickname,xtitle = textoidl('s = j/j_{tot}'),ytitle = 'P(s)';, charthick=2                                                                                                                                               
if n_elements(hs) lt n_elements(hg) then b = n_elements(hg) else b = n_elements(hs)
oplot, indgen(b)*0.05,(gmass/max(gmass))*(fdisk/.2121), linestyle=2
legend, ['DM', 'Gas'], linestyle=[0,2], /top, /right ;charthick = 2                                                                          
;xyouts, 1.8, -0.2, textoidl('s = j/j_{tot}'), charthick=3,
;charsize=1.5                                                                                                                
;xyouts, -.5, 0.45, 'P(s)', charthick=3, orientation=90, charsize=1.5                                                                                                                   
device, /close
end
