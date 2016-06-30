pro halo_sfh_master
dir = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/steps/h516.cosmo25cmb.3072g14HBWK.00080.dir'
file = 'h516.cosmo25cmb.3072g14HBWK.00080'
pfile = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK.param'
units = tipsyunits(pfile)
halo = [7,8,9,10,11,12,13,14,15,16,17,18]
halo_sfh,dir,file,halo,units,outplot = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g14HBWK/h516.cosmo25cmb.3072g14HBWK_small_halo_SFH.eps'
end


pro halo_sfh,dir,file,halo,units,outplot = outplot
!Y.STYLE = 1
!X.STYLE = 1
!P.THICK = 3.5
IF KEYWORD_SET(outplot) THEN BEGIN
    set_plot,'ps' 
    !P.CHARTHICK=4
    !X.THICK=4
    !Y.THICK=4
    !p.charsize=1.0
    !x.charsize=1.5;2.25
    !y.charsize=1.5;2.25
ENDIF ELSE BEGIN
    set_plot,'x'
    !P.CHARTHICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !p.charsize=1.0
    !x.charsize=1.5
    !y.charsize=1.5  
ENDELSE
!p.multi = 0

cd,dir
rtipsy,file,h,g,d,s
readarr,file + '.iord',h,iord,part = 'star',/ascii
readarr,file + '.amiga.grp',h,grp,part = 'star',/ascii

halo_ind = halo[0]
halo_stars_ind = where(grp eq halo_ind)
FOR i = 1,N_ELEMENTS(halo) - 1 DO BEGIN
    halo_ind = halo[i]
    halo_stars_ind = [halo_stars_ind, where(grp eq halo_ind)]
ENDFOR
halo_stars_iord = iord[halo_stars_ind]
halo_stars = s[halo_stars_ind]

;FOR i = 1, N_ELEMENTS(halo_stars_iord) DO BEGIN 
;END
IF (KEYWORD_SET(outplot)) THEN device,filename=outplot,/color,bits_per_pixel= 8,/times,xsize = 18,ysize = 12,xoffset =  2,yoffset = 2 ELSE window,0,xsize = 712,ysize = 39
sfr,halo_stars,massunits = units.massunit, timeunit = units.timeunit,xmargin = [14,3],ymargin = [8,3],xrange = [0,1],yrange = [0,0.04]
IF (KEYWORD_SET(outplot)) THEN device,/close
stop


end
