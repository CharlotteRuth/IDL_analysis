;Written for a homework assignment for PHY395, Spr 17 in which students create
;a CMD from SDSS data

PRO SDSS_CMD
loadct,39
H0 = 70
c = 3e5

;readcol,'CMD.txt',objid,mag_g,mag_i,mag_r,mag_u,z
;readcol,'CMD1000.txt',objid,mag_g,mag_i,mag_r,mag_u,z
readcol,'CMDmorph.txt',M_objid,M_mag_g,M_mag_i,M_mag_r,M_mag_u,M_z,ellip,spiral
M_i = mag_i - 5*alog10(z*c/H0*1e6/10)
M_r = mag_r - 5*alog10(z*c/H0*1e6/10)

levels = findgen(10)*500 + 500

;contour_plus,M_r,mag_g - mag_r,xmax = -19,ymax = 2,xmin = -22,ymin = 0.5,nlevels = 20,threshold = 250,xrange = [-19,-22],ytitle = 'g-r',xtitle = 'M_r'
contour_plus,M_r,mag_g - mag_r,xmax = -18,ymax = 2.5,xmin = -24,ymin = 0.0001,nlevels = 20,threshold = 250,xrange = [-18,-24],yrange = [0.2,1.2],ytitle = 'g-r',xtitle = 'M_r'
plot,M_r,mag_g -mag_r,psym = 3,xrange = [-18, -24],yrange = [0,2.5]
plot,M_r,mag_g -mag_r,psym = 3,xrange = [-18, -24],yrange = [0.2,1.5]

readcol,'CMDmorph.txt',M_objid,M_mag_g,M_mag_i,M_mag_r,M_mag_u,M_z,ellip,spiral
M_M_i = M_mag_i - 5*alog10(M_z*c/H0*1e6/10)
M_M_r = M_mag_r - 5*alog10(M_z*c/H0*1e6/10)
plot,M_M_r,M_mag_g -M_mag_r,psym = 3,xrange = [-18, -24],yrange = [0.2,1.6]
contour_plus,M_M_r,M_mag_g - M_mag_r,xmax = -18,xmin = -24,ymin = 0.2,ymax = 1.5,nlevels = 20,threshold = 250,xrange = [-18,-24],yrange = [0.2,1.5],ytitle = 'g-r',xtitle = 'M_r'
oplot,M_M_r,M_mag_g -M_mag_r,psym = 3,color = 100
oplot,M_M_r[where(spiral)],M_mag_g[where(spiral)] - M_mag_r[where(spiral)],psym = 3,color = 60
oplot,M_M_r[where(ellip)],M_mag_g[where(ellip)] - M_mag_r[where(ellip)],psym = 3,color = 254

levels = findgen(10)*250 + 250
plot,M_M_r,M_mag_u -M_mag_r,psym = 3,xrange = [-18, -24],yrange = [0.000001,5]
contour_plus,M_M_r,M_mag_u - M_mag_r,xmax = -20,xmin = -24,ymin = 0.5,ymax = 4.5,nlevels = 10,threshold = 250,xrange = [-20,-24],yrange = [0.5,4.5],ytitle = 'u-r',xtitle = 'M_r',export = hist2d
levels = findgen(10)*500 + 500
oplot,M_M_r[where(spiral)],M_mag_u[where(spiral)] - M_mag_r[where(spiral)],psym = 3,color = 60
oplot,M_M_r[where(ellip)],M_mag_u[where(ellip)] - M_mag_r[where(ellip)],psym = 3,color = 254

set_plot,'ps'
formatplot,outplot = 'CMD.eps'
device,filename = 'CMD.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2
contour_plus,M_M_r,M_mag_u - M_mag_r,xmax = -20,xmin = -24,ymin = 0.5,ymax = 4.5,nlevels = 10,threshold = 250,xrange = [-20,-24],yrange = [0.5,4.5],ytitle = 'u-r',xtitle = 'M_r',export = hist2d
oplot,M_M_r[where(spiral)],M_mag_u[where(spiral)] - M_mag_r[where(spiral)],psym = 3,color = 60
oplot,M_M_r[where(ellip)],M_mag_u[where(ellip)] - M_mag_r[where(ellip)],psym = 3,color = 254
device,/close
END

;----------------------------------
SELECT TOP 1000
    objid, dered_g, dered_i, dered_r, z as redshift   
    FROM SpecPhoto
    WHERE
    (class = 'GALAXY') AND
    z > 0.1 

SELECT
    g.objid, g.dered_g as g, g.dered_i as i, g.dered_r as r, g.dered_u as u, g.z as redshift,
    zns.elliptical,
    zns.spiral
FROM SpecPhoto as G
    JOIN ZooSpec AS zns
    ON g.objid = zns.objid
    WHERE (class = 'GALAXY') AND zns.nvote >= 10 AND g.z > 0.1
