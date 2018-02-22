PRO analyzeHIcubes,dir,cfile,thermal = thermal,physical_coord = physical_coord, angle = angle, dwarf = dwarf, doplot = doplot,xyrange = xyrange
telres_THINGS = 10.0 ;arcseconds
distance_THINGS = 5.0 ;Mpc
noise_THINGS = 0.65 ;miliJansky/beam, true for 1.3 and 2.6 km/s velocity resolution
IF keyword_set(dwarf) THEN vres_THINGS = 1.3 ELSE vres_THINGS = 2.6 ;km/s
limit_THINGS = 2.36e5*distance_THINGS^2*2.0*noise_THINGS/1000.0*vres_THINGS


IF NOT keyword_set(angle) THEN angle = 90
angle_str = strtrim(angle,2) + '.';+ 'KS.'

cd,dir
IF keyword_set(physical_coord) THEN BEGIN
    rtipsy,cfile,h,g,d,s,/justhead
;    header.CRVAL1 = header.CRVAL1;*h.time
;    header.CRVAL2 = header.CRVAL2;*h.time
;    header.CRVAL3 = header.CRVAL3;*h.time
;    header.CDELT1 = header.CDELT1;*h.time
;    header.CDELT2 = header.CDELT2;*h.time
;    header.CDELT3 = header.CDELT3;*h.time
    expansion = h.time;
    print,'Expansion: ',h.time
ENDIF

print,'*********************************************'
print,dir+'/'+cfile+'.cube.' + angle_str + 'fits'
print,'*********************************************'
cube = read_cube_fits(dir+'/'+cfile+'.cube.' + angle_str + 'fits',header,expansion = expansion)
shead = convertheader(header,max(cube),min(cube)) ;Converts header to strings in correct format to output to a fits file through MWRFITS
IF 1 THEN BEGIN
    IF angle NE 90 THEN BEGIN
        total_line_profile,cube,header,vaxis,spectrum,w,outfile = cfile,doplot = doplot,xyrange = xyrange ;Calculate line width as if from single dish data (i.e., don't smooth and remove sensetivity as if from THINGS)
        width = get_width_fit(vaxis, spectrum, /doplot, gal = cfile,incl=angle, pa=0) ;,sparange = sparange)
        print,'Width [km/s]:',width
        print,'V_max [km/s]:',width/2/sin(angle*!pi/180)

        openw,1,'linewidths.txt'
        printf,1,'Width [km/s]             Width/2 [km/s]             V_max [km/s]'
        printf,1,width,width/2,width/2/sin(angle*!pi/180)
        close,1
    ENDIF
ENDIF
moments,cube,header,m0,m1,m2,outplot = cfile + '.' + angle_str + 'cube',/fits;,/doplot ;Unsmoothed and sensetivity-limited moments

IF NOT (file_test(cfile+'.cube.' + angle_str + 'smoothed.fits')) THEN BEGIN
    smoothed = smooth_cube(cube,header,telres1 = telres_THINGS, telres2 = telres_THINGS, distance = distance_THINGS)
    cube = 0
    smoothed[where(smoothed lt limit_THINGS)] = 0 ;take into account resolution
    smoothed = smoothed*1.0/(1.133*10*10)*(1.5)^2; change back to total solar masses rather than solar masses per beam
    shead = convertheader(header,max(smoothed),min(smoothed))
    mwrfits,smoothed,cfile+'.cube.' + angle_str + 'smoothed.fits',shead
ENDIF ELSE BEGIN
    cube = 0    
    smoothed = read_cube_fits(cfile+'.cube.' + angle_str + 'smoothed.fits',header)
ENDELSE
zeropoint = MIN(smoothed)
smoothed[where(smoothed eq zeropoint)] = 0
;look for velocity dispersions around 8-10 km/s
moments,smoothed,header,m0,m1,m2,outplot = cfile + '.' + angle_str + 'smoothed',/fits;,/doplot ;THINGS-like moments
;total_line_profile,smoothed,header,vaxis,spectrum,w,outfile = file
smoothed = 0
END

PRO analyzeHIcubes_master

base = '/astro/net/scratch2/christensen/MolecH/Cosmo/'
;directories = ['h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir', $
;directories = ['h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00168.dir', $
directories = ['h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00512.dir', $
               'h516.cosmo25cmb.1536g2HBWK/steps/h516.cosmo25cmb.1536g2HBWK.00168.dir', $
               'h516.cosmo25cmb.1536g2HBWK_nocool/steps/h516.cosmo25cmb.1536g2HBWK.00168.dir']
base = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/'
directories = ['h516.cosmo25cmb.1536g1MBWK/steps/h516.cosmo25cmb.1536g1MBWK.00512.dir', $
               'h516.cosmo25cmb.1536g14HBWK/steps/h516.cosmo25cmb.1536g14HBWK.00512.dir']
base = '/astro/store/student-scratch1/christensen/MolecH/Cosmo/'
directories = ['h516.cosmo25cmb.2304g/h516.cosmo25cmb.2304g14HBWK/steps/h516.cosmo25cmb.2304g14HBWK.00512.dir',$
               'h516.cosmo25cmb.3072g/h516.cosmo25cmb.3072g1MBWK/steps/h516.cosmo25cmb.3072g1MBWK.00492.dir']

base = '/home/byrnelin/MAP/'
directories = ['Disk_Collapse_1e6_H2_JC','Disk_Collapse_1e6_TC_JC','Disk_Collapse_1e6_JC','Disk_Collapse_1e6_LW2']

directories = base + directories
;files = ['h516.cosmo25cmb.1536g1MBWK.00512', $
;files = ['h516.cosmo25cmb.1536g1MBWK.00168', $
files = ['h516.cosmo25cmb.1536g2HBWK.00512', $
         'h516.cosmo25cmb.1536g2HBWK.00168', $
         'h516.cosmo25cmb.1536g2HBWK.00168']
files = ['h516.cosmo25cmb.1536g1MBWK.00512', $
         'h516.cosmo25cmb.1536g14HBWK.00512']
files = ['h516.cosmo25cmb.2304g14HBWK.00512',$
         'h516.cosmo25cmb.3072g1MBWK.00492']
files = files + '.halo.1.cube'

files = ['Disk_Collapse_1e6_H2.00100.scalez1test.000100','Disk_Collapse_1e6.00100.scalez1test.000100','Disk_Collapse_1e6_H2_Shielding.00100.scalez1test.000100','Disk_Collapse_1e6.00100.scalez1test.000100']

FOR i = 0, N_ELEMENTS(files) - 1 DO BEGIN
   analyzeHIcubes,directories[i],files[i],angle = 45
   analyzeHIcubes,directories[i],files[i],angle = 90
ENDFOR
END
