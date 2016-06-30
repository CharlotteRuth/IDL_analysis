;plot the profile of the galaxy, with the galfit fits overplotted
;input the galfit fits 'galfit.0x' 
;and fov,npix of your sunrise image, as set in mcrx.stub
pro checkfit,galfit,fov,npix

prefix='./'
file = prefix+'gal.fits'
fitfile=prefix+galfit
readcol,fitfile,dummy,fits,/silent,format='A,F'
readcol,fitfile,dummy,com1,com2,/silent,format='A,F,F'
zeropoint=fits[3]


;SB profile of the image, using circular annuli
gal =mrdfits(file)
makex, gal,x,y  ;set up coordinate grid
pixKpc = 1.0*fov/(1.0*npix)
sbfactor=10^(-zeropoint/2.5)
cx=(com1[3]+com1[12]-npix)/2
cy=(com2[3]+com2[12]-npix)/2

radius = sqrt((x)^2+(y)^2) ; total face on
hist = hist1d(radius,gal,obin=xvals,binsize=2) 
numh = hist1d(radius,binsize=2)
SBprofile = hist/numh
rad_kpc = xvals*pixKpc
profile = -2.5*alog10(sbfactor*SBprofile)
;****************************************************
; galfit fits 
k = 1.9992*fits[10] - 0.3271

ftot=10^((-fits[8]+zeropoint)/2.5)
Rb=!pi*(fits[15]+2)/(4*BETA(1/(fits[15]+2),1+1/(fits[15]+2)))
ue=ftot*Rb/(2*!pi*fits[9]^2*exp(k)*k^(-2*fits[10])*GAMMA(2*fits[10])*fits[13])*sbfactor
bulgefit = abs(ue)*exp(-k*((rad_kpc/fits[9]/pixKpc)^(1/fits[10])-1))
bulgeplot = -2.5*alog10(bulgefit)

ftot=10^((-fits[18]+zeropoint)/2.5)
Rb=!pi*(fits[25]+2)/(4*BETA(1/(fits[25]+2),1+1/(fits[25]+2)))
uo=ftot*Rb/(2*!pi*fits[19]^2*fits[23])*sbfactor
diskfit = abs(uo)*exp(-(rad_kpc/fits[19]/pixKpc))
diskplot = -2.5*alog10(diskfit)
fit=bulgefit+diskfit
fitplot = -2.5*alog10(fit)

;************************************************
;plotting
set_plot,'x'
set_plot,'ps'
device,/portrait,filename='./checkfit.eps', /color

xmin=0.
xmax=10.
region=where(rad_kpc gt xmin and rad_kpc lt xmax)
ymin=max(profile[region])+0.8
ymax=min(profile[region])-0.8
ymin = 33
pixmin=xmin*(npix/fov)
pixmax=xmax*(npix/fov)


loadct,39
plot,rad_kpc,profile,xrange=[xmin,xmax],/xstyle,thick=2,charsize=1.5,xtit='Radius (kpc)',yrange=[ymin,ymax],/ystyle,psym=0, ticklen=0;,title = 'Scale Length = '+STRTRIM(fits[19])+' kpc'
AXIS, XAXIS=1, XRANGE=[pixmin,pixmax],  xticks=1, xtickinterval=50, ticklen=.02, /xstyle,charsize=1.5
AXIS, XAXIS=0,xrange=[xmin,xmax],/xstyle,charsize=1.5
AXIS, YAXIS=0,yrange=[ymin,ymax],/ystyle,charsize=1.5,ytitle='Magnitude'
;yrange=[22.49,16.01]

oplot,rad_kpc,diskplot,linestyle=1,color=150
oplot,rad_kpc,bulgeplot,linestyle=2,color=250
oplot,rad_kpc,fitplot,linestyle=3,color=50
legend, ['Sersic', 'Expdisk', 'Sersic+Expdisk','B/D = '+strtrim(10^(-fits[8]/2.5)/10^(-fits[18]/2.5),2)], linestyle=[2,1,3,-0],charsize=1.5,color=[250,150,50,0],/right
device,/close

set_plot,'x'
window,1
plot,rad_kpc,profile,xrange=[xmin,xmax],/xstyle,thick=2,charsize=1.5,xtit='Radius (kpc)',yrange=[ymin,ymax],/ystyle,psym=0, ticklen=0;,title = 'Scale Length = '+STRTRIM(fits[19])+' kpc'
AXIS, XAXIS=1, XRANGE=[pixmin,pixmax],  xticks=1, xtickinterval=50, ticklen=.02, /xstyle,charsize=1.5
AXIS, XAXIS=0,xrange=[xmin,xmax],/xstyle,charsize=1.5
AXIS, YAXIS=0,yrange=[ymin,ymax],/ystyle,charsize=1.5,ytitle='Magnitude'
;yrange=[22.49,16.01]

oplot,rad_kpc,diskplot,linestyle=1,color=150
oplot,rad_kpc,bulgeplot,linestyle=2,color=250
oplot,rad_kpc,fitplot,linestyle=3,color=50
legend, ['Sersic', 'Expdisk', 'Sersic+Expdisk','B/D = '+strtrim(10^(-fits[8]/2.5)/10^(-fits[18]/2.5),2)], linestyle=[2,1,3,-0],charsize=1.5,color=[250,150,50,0],/right
stop

print,'scalelength =',fits[19],' kpc'
print,'B/D =',10^(-fits[8]/2.5)/10^(-fits[18]/2.5)
stop
end

