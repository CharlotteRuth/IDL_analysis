PRO plots

a=findgen(17)*(!PI/8.)
usersym,cos(a),sin(a),/fill


readcol, "galaxy_info.tbl",format='a,f,f,f,f,f,f,f',bandpass, lambda,galfit, disk_scalelength,BD, mu_0, bulge_n, bulge_mag, disk_mag

readcol, "galaxy_info_3comp.tbl",format='a,f,f,f,f,f,f,f,f,f,f',bandpass3, lambda3, galfit3, disk_scalelength3, BD3, mu_03, bulge_n3, bar_n3, bulge_mag3, disk_3mag, bar_mag

;set_plot,'ps'
;device,portrait=1,filename='galfit_F606.eps',/color,/encapsulate


x=lambda
y=BD
y2=BD3

loadct,12
plot,x,y,xrange=[min(x),max(x)],/xstyle,/xlog,thick=2,charsize=1.3,xtit='Wavelength (nm)',yrange=[min(y),min(y)],/ystyle,/nodata, ticklen=0, pos=[.2,.2,.9,.85]

AXIS, YAXIS=0,yrange=[ymin,ymax],/ystyle,charsize=1.3,ytitle='Bulge/Disk Ratio'

oplot, y, psym=6,color=50,symsize=1.1,thick=3
oplot, y2,psym=8,color=250,symsize=.85


;legend, ['Bulge', 'Disk', 'Bar', 'Total'], linestyle=[2,1,0,3],color=[250,150,200,100], charsize=1.5, /top, /right



;device,/close
;set_plot,'x'


return
end
