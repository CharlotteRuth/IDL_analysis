pro coldparticles
loadct,39
cm_in_kpc = 3.08568025e21

spawn,'ls *.param',pfilelist
pfile = pfilelist[0]
units = tipsyunits(pfile)

filebase = 'h516.cosmo25cmb.1536g8HBWK.JeansSF.00408'
steps = ['00001','00002','00003']
steps = ['00002','00003']
indcold = [192039,210846,411405,607324,661804,715958]
colors = (FINDGEN(N_ELEMENTS(indcold))+1)/N_ELEMENTS(indcold)*240
FOR i=0,N_ELEMENTS(steps) - 1 DO BEGIN
    name = filebase + '.' + steps[i]
    window,i
    !p.multi = [0,2,1]
    rtipsy,name,h,g,d,s
    plot,g.dens*units.rhounit,g.tempg,psym = 3,/xlog,/ylog,xrange = [1e-10,1e2],yrange = [1,1e8]
    FOR ip=0,N_ELEMENTS(indcold) - 1 DO oplot,[g[indcold[ip]].dens,g[indcold[ip]].dens]*units.rhounit,[g[indcold[ip]].tempg,g[indcold[ip]].tempg],psym = 2,color = colors[ip]
    indden = where(g.dens*units.rhounit gt 1e-8)
    ghighden = g[indden]
    highden = alog10(g[indden].dens*units.rhounit)
    dencolors = (highden+8)/10.0*240
    loadct,0
    plot,g[indden].x,g[indden].y,psym = 3,xrange = [-0.04,-0.01],yrange = [0.05,0.07];first step, all particles
;    plot,g[indden].x,g[indden].y,psym = 3,xrange = [-0.020,-0.015],yrange = [0.05,0.06] ;third step, all particles
    plot,g[indden].x,g[indden].y,psym = 3,xrange = [-0.01555,-0.0154],yrange = [0.0518,0.05195];third step, disk
;    plot,g[indden].x,g[indden].y,psym = 3,xrange = [-0.01555,-0.0154],yrange = [0.0518,0.0520]
;    plot,g[indden].x,g[indden].y,psym = 3,xrange = [-0.01551,-0.01545],yrange = [0.05180,0.05182]
    FOR ip = 0L,N_ELEMENTS(highden)-1 DO oplot,[ghighden[ip].x,ghighden[ip].x],[ghighden[ip].y,ghighden[ip].y],psym = 3,color = dencolors[ip]
    loadct,39
    FOR ip=0,N_ELEMENTS(indcold) - 1 DO oplot,[g[indcold[ip]].x,g[indcold[ip]].x],[g[indcold[ip]].y,g[indcold[ip]].y],psym = 2,color = colors[ip] 
    stop
    !p.multi = 0
ENDFOR

stop
end
