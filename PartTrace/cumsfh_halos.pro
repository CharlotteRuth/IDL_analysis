PRO MASTER
;dir='/nobackupp8/crchrist/MolecH/cptmarvel.cosmo25.4096g/cptmarvel.cosmo25cmb.4096g5HbwK1BH/'
;filebase='cptmarvel.cosmo25cmb.4096g5HbwK1BH'
;finalids = ['1','2','4','5','6','7','10','11','13','14']

;dir='/nobackupp8/crchrist/MolecH/rogue.cosmo25cmb.4096g/rogue.cosmo25cmb.4096g5HbwK1BH/'
;filebase='rogue.cosmo25cmb.4096g5HbwK1BH'
;finalids = ['1','3','7','8','10','11','12','16','17','18','30','31','32','34']

;dir = '/nobackupp8/crchrist/MolecH/elektra.cosmo25cmb.4096g/elektra.cosmo25cmb.4096g5HbwK1BH/'
;filebase='elektra.cosmo25cmb.4096g5HbwK1BH'
;finalids=['1','2','3','4','5','8','9','10','11','12','17','18'];8,10,12

dir = '/nobackupp8/crchrist/MolecH/storm.cosmo25cmb.4096g/storm.cosmo25cmb.4096g5HbwK1BH/' 
filebase='storm.cosmo25cmb.4096g5HbwK1BH' 
finalids= ['1','2','3','4','5','6','7','8','10','11','13','14','15','16','17','23','24','28','34','35','43','49','50','60'] 

laststep = '004096'
cumsfh_halos,dir,filebase,finalids, laststep,/outplot

END

PRO cumsfh_halos, dir, filebase, finalids, laststep,outplot = outplot
formatplot,outplot = outplot
loadct,39
colors = findgen(n_elements(finalids))/n_elements(finalids)*(256 - 13) + 13
IF keyword_set(outplot) THEN BEGIN
    xsize = 18
    ysize = 12
ENDIF ELSE BEGIN
    xsize = 800
    ysize = 500
ENDELSE

cd,dir
rtipsy,filebase+'.'+laststep+'/'+filebase+'.'+laststep,h,g,d,s,/justhead
units = tipsyunits(filebase + '.param')
grp = read_lon_array(filebase+'.'+laststep+'/'+filebase+'.'+laststep+'.amiga.grp')
grpstar = grp[h.ngas+h.ndark:n_elements(grp) - 1]
iord = read_lon_array(filebase+'.'+laststep+'/'+filebase+'.'+laststep+'.iord')
iordstar = iord[h.ngas+h.ndark:n_elements(iord) - 1]
sl = rstarlog(filebase+'.starlog',/molecularBH);,/BIG)
match2,sl.iorderstar,iordstar,ind1,ind2
IF n_elements(where(ind1 NE -1)) LT 10 THEN BEGIN
    sl = rstarlog(filebase+'.starlog',/molecularH,/BIG)
    match2,sl.iorderstar,iordstar,ind1,ind2
ENDIF
sltrim = sl[where(ind1 NE -1)]

IF keyword_set(outplot) THEN device,filename = filebase + '.cumsfh.eps',/encapsulated,/color,/times,ysize=ysize,xsize=xsize,bits_per_pixel= 8 ELSE window,0,ysize=ysize,xsize=xsize
FOR ih = 0, n_elements(finalids) - 1 DO BEGIN
    slhalo = sltrim[where(grpstar EQ finalids[ih])]
    slhalo = slhalo[where(slhalo.timeform GT 0)]
    IF ih EQ 0 THEN histogramp,slhalo.timeform*units.timeunit/1e9,weight = slhalo.massform*units.massunit,/normalize,/cum,min = 0,max = 13.8,nbins = 200,title = filebase,xtitle = 'Time [Gyr]',ytitle = 'Cumulative SFH',yrange = [0,1]
    loadcv,80,/reverse
    histogramp,slhalo.timeform*units.timeunit/1e9,weight = slhalo.massform*units.massunit,/normalize,/cum,min = 0,max = 13.8,nbins = 200,/overplot,color = colors[ih],thick = 3
ENDFOR
IF keyword_set(outplot) THEN device,/close ELSE stop
END
