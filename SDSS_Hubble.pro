;Written for a homework assignment for PHY395, Spr 17 in which students create
;a CMD from SDSS data
readcol,'HubbleMorph.txt',M_objid,M_mag_g,M_mag_i,M_mag_r,M_mag_u,M_z,ellip,spiral

Mr0 = -21.5
distance = 10*10^((M_mag_r[where(spiral)] - Mr0)/5)/1e6 ;Mpc
vel = M_z[where(spiral)]*c ;km/s

measure_errors = (findgen(n_elements(vel)) + 1)*1e4
result = linfit([0,distance], [0,vel], MEASURE_ERRORS=[1e-6,measure_errors])

nlevels = 10
set_plot,'ps'
formatplot,outplot = 'SDSSHubble.eps'
device,filename = 'SDSSHubble.eps',/color,bits_per_pixel= 8,xsize = xsize,ysize = ysize,xoffset =  2,yoffset =  2
contour_plus,distance,vel,xmax = 800,ymax = 8e4,xmin = 0.0001,ymin = 0.0001,nlevels = nlevels,threshold = 20,xrange = [0,800],yrange = [0,8e4],ytitle = 'Velocity [km/s]',xtitle = 'Distance [Mpc]',xbinsize = 10,ybinsize = 1000,C_COLOR = (findgen(nlevels) + 1)*254/nlevels,title = textoidl('H_0 = '+strtrim(string(round(result[1])),1)+' km/s/Mpc')
oplot,[0,800],[result[0],result[0] + 800*result[1]]
device,/close

SELECT
    g.objid, g.dered_g as g, g.dered_i as i, g.dered_r as r, g.dered_u as u, g.z as redshift,
    zns.elliptical,
    zns.spiral
FROM SpecPhoto as G
    JOIN ZooSpec AS zns
    ON g.objid = zns.objid
    WHERE (class = 'GALAXY') AND zns.nvote >= 10 AND zns.spiral = 1
