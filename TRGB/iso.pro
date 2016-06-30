;extrap_iso is in extrap_trgb.pro

function findtip,xvalues,yvalues
;Find the tip of the RGB
previousx = xvalues[0]
previousy = yvalues[0]
for i = 1, (size(xvalues))[2] - 1 DO BEGIN
    if (previousx gt xvalues[i]) AND (previousy lt yvalues[i]-0.5) then return,i-1
    previousx = xvalues[i]
    previousy = yvalues[i]
endfor
return,i
end

function readin,filename
;This reads in the isochrones and then returns the V and I magnitudes
;for the isochrone upto the tip of the RGB
filename = '/astro/net/scratch2/christensen/TRGBDust/isochrones/customized/isoc_z0' + filename + '.dat'
ncol = 16 
nrow = 184
header = strarr(2)
data = dblarr(ncol,nrow)

close,1
openr,1,filename
readf,1,header
readf,1,data
close,1
mags = dblarr(2,nrow)
mags[0,*] = data[9,*] ;V magnitude
mags[1,*] = data[11,*] ;I magnitude
tip = findtip(mags[0,*]-mags[1,*],mags[1,*])
mags = EXTRAC(mags,0,0,2,tip)
return, mags
end

FUNCTION findcolor,mag,tips
delta = 1
FOR i=0,N_ELEMENTS(tips[2,*])-1 DO BEGIN
    IF ABS(tips[2,i] - mag) LT delta THEN BEGIN
        color = tips[1,i]
        delta = ABS(tips[2,i] - mag)
;        print,i,": ",delta
    ENDIF
ENDFOR
RETURN,color
END

pro iso
;Charlotte Christensen, 7/17/07
;The program reads in isochrones of different metallicities
;addx = -0.4
obsdata_file = '/astro/net/scratch2/christensen/TRGBDust/NGC0253-WIDE1/'
iso_file = '/astro/net/scratch2/christensen/TRGBDust/isochrones/metalgrid/'


;;;;;;;;;;;;;;;;;;;;;;;;READ IN OBSERVATIONAL DATA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
addy = 26.8+0.88
galaxy = 'NGC0253-WIDE1'
infile=obsdata_file + galaxy + '.vi2.fits'
columns = ['mag1_acs', 'mag1_err','mag2_acs','mag2_err']
data = mrdfits(infile,1,columns=columns)
mI = data.mag2_acs ;i
mColor = data.mag1_acs - data.mag2_acs ;v - i
loadct,39
set_plot,'x'
;set_plot,'ps'

;;;;;;;;;;;;;;;;;;;;;;READ IN ISOCHRONES;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
filelength = 120 ;Find using wc on it (300)or decide what range is desired
num_isos = 12
isos_name_init = strarr(filelength)
infile_iso = iso_file + 'index.txt'
close,1
openr,1,infile_iso
readf,1,isos_name_init
close,1
isos_names = Reform((Reform(isos_name_init,filelength/num_isos,num_isos))[0,*])
elements = 150
isochrones = dblarr(3,elements,num_isos);[quality, points, isochrones]


FOR i = 0, (num_isos - 1) DO BEGIN
    file = extrap_iso(iso_file + isos_names[i],elements)
    isochrones[0,*,i] = file.z                      ;metalicity
    isochrones[1,*,i] = file.f606mag - file.f814mag ;color
    isochrones[2,*,i] = file.f814mag + addy                ;imagnitude
ENDFOR

READCOL,iso_file+"Tips.dat",z,f475mag,f606mag,f814mag,FORMAT='F,X,X,X,X,X,X,X,X,F,X,X,F,X,X,X,X,F',/SILENT
tips = dblarr(3,filelength)
tips[0,*] = z[0:filelength-1]
tips[1,*] = f606mag[0:filelength-1] - f814mag[0:filelength-1]
tips[2,*] = f814mag[0:filelength-1]+addy

;;;;;;;;;;;;;;; Finding the tip for each of the segments;;;;;;;;;;;;;;;;;;;;
mI = data.mag2_acs
mColor = data.mag1_acs - data.mag2_acs
length = N_ELEMENTS(mI)
trgb_result = fltarr(num_isos-1,6) ; Structures to hold results
openw,1,'iso_tips.dat'
printf,1,'First Isochrone ','Second Isochrone ','TRGB Color','Error in Color ','Magnitude ','Error in Magnitude',' Ave Magnitude of Tip Points'
FOR ct = 0, num_isos-2 DO BEGIN
    print,"Metalicity: ",isos_names[ct]
    length = N_ELEMENTS(mI)
    between = replicate(0,length)
    FOR i = 0L, length - 1 DO BEGIN
        IF ((below(isochrones[1,*,ct],isochrones[2,*,ct],mColor[i],mI[i]) EQ 1) AND (above(isochrones[1,*,ct+1],isochrones[2,*,ct+1],mColor[i],mI[i]) EQ 1)) THEN between[i] = 1 
;Determins whether a star is between two isochrones
    ENDFOR
    between = Where(between EQ 1)
    window,1
    !P.MULTI=0
    plot,mColor,mI,xtitle = 'V - I', ytitle = 'I', title = 'Isochrone CMDs',yrange=[25.5,22],xrange=[-0.5,4.0],psym = 3
    FOR i = 0, (num_isos - 1) DO BEGIN
        oplot,isochrones[1,*,i],isochrones[2,*,i],color = 220 - i*(220./num_isos)
    ENDFOR
    oplot,mColor[between], mI[between],color = 230,psym = 3

    result = trgb(galaxy,(data.mag1_acs)[between],(data.mag1_err)[between],(data.mag2_acs)[between],(data.mag2_err)[between],tips[1:2,ct],tips[1:2,ct+1],isochrones[1:2,*,ct],isochrones[1:2,*,ct+1],ct,/usevi)

    printf,1,isos_names[ct]," ",isos_names[ct+1],result[0],result[1],result[2],result[3]
ENDFOR
close,1
END
