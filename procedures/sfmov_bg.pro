pro sfmov_bg,tipsyfile,ppmroot
;makes a movie of star formation with background galaxy.

windowsize=500
window,xsize=windowsize,ysize=windowsize
paramfile = strmid(tipsyfile,0,strlen(tipsyfile)-5) + "param"
dKpcUnit = 1e5
a = timeunits(paramfile)
dtUnit = a[0]
dt  = 985624.
twidth = 1e7
rformfile=tipsyfile+".rform"
rdfloat,rformfile,rform,skip=1
rtipsy,tipsyfile,h,g,d,s
ns = h.nstar
ntot = h.n
ngd= h.ngas + h.ndark

tform = s.tform * dtUnit

ind = where(tform gt 0)
mint = min(tform[ind])
maxt = max(tform[ind])

r = replicate({x:0.,y:0.,z:0.},ns)
r.x = -rform[ngd:ntot-1]
r.y = rform[ntot+ngd:2*ntot-1]
r.z = rform[2*ntot+ngd:3*ntot-1]

!ORDER=1
 ;mpegid = mpeg_open([400,400]) 
 
nframes=FIX((maxt-mint)/dt)
print,nframes
outim=INDGEN(3,windowsize,windowsize)

for i=0,nframes-1 DO BEGIN
     t=mint+FLOAT(i)*dt
     

     ;set up plotting window
     plot,FINDGEN(20),xrange=[-20,20],yrange=[-20,20],/nodata,/isotropic,xtit='X [kpc]',ytit = 'Y [kpc]',xthick=2,ythick=2

     ;display the ppm image
     ppmext='000000000'+STRTRIM(i+1,2)
     ppmin=ppmroot+"."+STRMID(ppmext,8,/reverse_offset)+".ppm"
     read_ppm,ppmin,inimage
     my_imgunder,REFORM(inimage[0,*,*]),channel=1
     my_imgunder,REFORM(inimage[1,*,*]),channel=2
     my_imgunder,REFORM(inimage[2,*,*]),channel=3

     labelt=FIX(FIX(i/10.)*10.)
     label='t='+STRTRIM(labelt,2)+' Myr'
     xyouts,5,15,label,charsize=2.0,charthick=2.0

     in = where( (tform gt t-twidth) AND (tform lt t+twidth),nin )
     if (nin gt 0) then begin         
         plot, r[in].x*dKpcUnit, r[in].y*dKpcUnit,psym=4,/isotropic,xtit='X [kpc]',ytit = 'Y [kpc]',xrange=[-20,20],yrange=[-20,20],/noerase,/nodata,xthick=2,ythick=2
         loadct,2,/silent
         oplot, r[in].x*dKpcUnit, r[in].y*dKpcUnit,psym=4,color=25,symsize=0.5
         loadct,0,/silent
     ENDIF

     if (i LT 10) then begin
         ppmout = "sf-x000"+strcompress(string(i),/remove_all)+".ppm"
     endif else if (i LT 100) then begin
         ppmout = "sf-x00"+strcompress(string(i),/remove_all)+".ppm"
     endif else if (i LT 1000) then begin
         ppmout = "sf-x0"+strcompress(string(i),/remove_all)+".ppm"
     endif else ppmout = "sf-x"+strcompress(string(i),/remove_all)+".ppm"
     
     outim[0,*,*] = tvrd(channel=1)
     outim[1,*,*] = tvrd(channel=2)
     outim[2,*,*] = tvrd(channel=3)
     
     write_ppm,ppmout,outim
         
                                ;     mpeg_put,mpegid,window=0,frame=fix(t/dt + 0.1)
endfor
!ORDER=0  
 ;mpeg_save,mpegid,filename='starform.mpg'
 ;mpeg_close,mpegid
end
