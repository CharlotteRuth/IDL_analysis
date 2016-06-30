;cd,'/astro/net/scratch2/christensen/MolecH/Cosmo/h516.cosmo25cmb.1536g/h516.cosmo25cmb.1536g3HBWK/steps/h516.cosmo25cmb.1536g3HBWK.00512.dir'
;writehimacro,'h516.cosmo25cmb.1536g3HBWK.00512.halo.1',pfile = '../../h516.cosmo25cmb.1536g3HBWK.param',/dwarfv
;writehimacro,'h516.cosmo25cmb.3072g14HBWK.00512.halo.2',pfile= '../../h516.cosmo25cmb.3072g14HBWK.param',/dwarfv,/molecularH

;---------- For THERMAL
;openb h516.cosmo25cmb.3072g14HBWK.00512.halo.1.std
;loads 1
;xall
;redshift   -3.7747583e-15 c 73 0.24 0.76 p 1
;coolconstants 25000 2.31000e+15 1
;setsphere 1 0 0 0 0.00047999999
;abox 1
;viewgas hfrac 0 1
;rotate up 45
;viewgas hfrac 0 1
;cube gas h516.cosmo25cmb.3072g14HBWK.00512.halo.1.thermal.cube.fits
;1024 -400 400 5.2 thermal

;dist_unit = 50000
;mass_unit = 1.84795e16

;dist_unit = 25000
;mass_unit = 2.310e15

;radius = 24
;writeHImacro, infile, dist_units = dist_units, mass_unit = mass_unit, pfile = pfile, radius = radius

;absv = 249.6

;This function writes a macro to make an HI data cube through tipsy
pro writehimacro, infile, dist_unit = dist_unit, mass_unit = mass_unit, pfile = pfile, radius = radius, molecularH = molecularH,SD_obs = SD_obs, dwarfv = dwarfv, npixels = npixels,fileext = fileext
;Use keyword 'dwarfv' if V_max is less than or equal to 83 km/s in
;physical units

IF keyword_set(dwarfv) THEN BEGIN
    absv = 83.2
    deltav = 1.3
ENDIF ELSE BEGIN
    absv = 166.4
    deltav = 2.6
ENDELSE
absv = 249.6

;spawn,'ls *.AHF_halos', z_file
;z_sec = strsplit(z_file,'.',/extract)
;z = strmid(z_sec[N_ELEMENTS(z_sec) - 3],1) + '.'+ z_sec[N_ELEMENTS(z_sec) - 2]
;time = z+1
rtipsy,infile,h,/justhead
z = 1/h.time - 1
absv_comove = strtrim(absv/h.time,2)
deltav_comove = strtrim(deltav/h.time,2)
IF keyword_set(pfile) THEN BEGIN
    unit = tipsyunits(pfile,silent = 1)
    dist_unit = unit.lengthunit
    mass_unit = unit.massunit
ENDIF
IF NOT keyword_set(radius) THEN radius = 18.62;kpc    used to be 12.0 before 7/26/11
IF NOT keyword_set(npixels) THEN npixels = 1024
IF NOT keyword_set(fileext) THEN fileext = ''

;close,/all
openw,1,infile+'.HI' + fileext + '.macro'
printf,1,"makeHIcube"
printf,1,""
IF keyword_set(SD_obs) THEN printf,1,"openb "+infile ELSE printf,1,"openb "+infile+".std"
printf,1,"loads 1"
printf,1,"xall"
printf,1,"redshift ",z," c 73 0.24 0.76 p 1"
printf,1,"coolconstants "+strtrim(dist_unit,2) + " " +strtrim(mass_unit,2)+" 1"
printf,1,"setsphere 1 0 0 0 "+strtrim(radius/dist_unit/h.time,2)
printf,1,"abox 1"
printf,1,"loadhneutral ascii "+infile+".HI"
printf,1,"viewgas hfrac 0 1"

IF keyword_set(SD_obs) THEN BEGIN
;Useful for creating HI-H2 transition plots
    printf,1,"rotate up 90"
    printf,1,"viewgas hfrac 0 1"
    printf,1,"cube gas "+infile+".cubeHI.fits 1024 -200 200 10"
    printf,1,"cube gas "+infile+".cubeH2.fits 1024 -200 200 10 molec"    
;    printf,1,"loadsparray ascii "+infile+".H2"
;    printf,1,"cube_array "+infile+".H2.fo.fits t log h 1024"
;    printf,1,"loadsparray ascii "+infile+".HI"
;    printf,1,"cube_array "+infile+".HI.fo.fits t log h 1024"
ENDIF ELSE BEGIN
    printf,1,"rotate up 45"
    printf,1,"viewgas hfrac 0 1"
    printf,1,"cube gas "+infile+".cube.45" + fileext + ".fits " + strtrim(npixels) + " -" + absv_comove + " " + absv_comove + " " + deltav_comove
    printf,1,"loadsparray ascii "+infile+".HI"
    printf,1,"cube_array "+infile+".HI.arr.45" + fileext + ".fits t log h 1024"
    IF keyword_set(molecularH) THEN  BEGIN
        printf,1,"loadsparray ascii "+infile+".H2"
        printf,1,"cube_array "+infile+".H2.arr.45" + fileext + ".fits t log h 1024"
    ENDIF
    printf,1,"zgas hfrac 0 1"
    printf,1,"viewgas hfrac 0 1"
    printf,1,"cube gas "+infile+".cube.90" + fileext + ".fits " + strtrim(npixels) + " -" + absv_comove + " " + absv_comove + " " + deltav_comove
    printf,1,"loadsparray ascii "+infile+".HI"
    printf,1,"cube_array "+infile+".HI.arr.90" + fileext + ".fits t log h 1024"
    IF keyword_set(molecularH) THEN  BEGIN
        printf,1,"loadsparray ascii "+infile+".H2"
        printf,1,"cube_array "+infile+".H2.arr.90" + fileext + ".fits t log h 1024"
    ENDIF
ENDELSE
printf,1,"closeb"
printf,1,"end"
close,1


spawn,'pwd',dir
;--------------------------------------------------

openw,1,'analyzeHIcubes.cfg'
printf,1,'getenv = true'
printf,1,'Universe = vanilla'
printf,1,'Executable = ' + dir + '/con_script.analyzeHIcubes'
printf,1,'Initialdir = ' + dir
printf,1,'output =    analyzeHIcubes.out'
printf,1,'error =     analyzeHIcubes.error'
printf,1,'Log =       analyzeHIcubes.log'
printf,1,' '
printf,1,'# condor requirements'
printf,1,'# this makes sure the job goes to elektra, and ensures that elektra'
printf,1,'# acknowledges the job as a valid SMP job'
printf,1,'+JobClass = "SMP"'
printf,1,'Requirements = MachineClass == "SMP"'
printf,1,'+CPUS = 1'
printf,1,' '
printf,1,'Queue'
close,1

;--------------------------------------------------

openw,1,'con_script.analyzeHIcubes'
printf,1,'#!/bin/bash'

printf,1,'export PATH=/astro/users/christensen/code/kcorrect/bin:/astro/condor/RH5.x86_64/bin:/astro/condor/RH5.x86_64/sbin:/astro/users/christensen/bin:/astro/net/jdk/bin:/astro/net/mozilla:/local/bin:/astro/net/python/bin/:/astro/apps/pkg/python/bin:/astro/net/intel/compiler91/bin:/astro/users/christensen/globus/bin:/astro/users/christensen/globus/sbin:/bin:/usr/bin:/usr/X11R6/bin:/usr/local/bin:/usr/sbin:/astro/apps:/astro/net/bin:/astro/net/condor/bin:/astro/users/christensen/Applications:.:/astro/net/condor/bin:/users/trq/bin.i386:/astro/users/patrik/bin'
printf,1,'export IDL_DIR=/astro/net/idl'
printf,1,'export LIB=/net/idllib'
printf,1,'export IDL_STARTUP=/astro/users/christensen/idl/idl_startup.pro'
printf,1,'export IDL_PATH=/astro/apps/idlpkg/sdssidl/pro:+/astro/users/christensen/code/IDL:+/usr/local/opt/itt/idl/lib:+/usr/local/opt/itt/idl/examples:+/astro/net/idllib:+/astro/users/christensen/code/idlutils/goddard/pro:+/astro/users/christensen/code/idlutils/pro:+/astro/apps/idlpkg/sdssidl/pro:+/astro/users/christensen/code/kcorrect/pro'
str_dir = "'" + dir + "'"
printf,1,'echo "dir = ' +  str_dir + '" > analyzeHIcubes.idl.batch'
str_file = "'" + infile + "'"
printf,1,'echo "file = ' + str_file + '" >> analyzeHIcubes.idl.batch'
printf,1,'echo ".r /astro/users/christensen/code/IDL/HIcubes/convertheader.pro" >> analyzeHIcubes.idl.batch'
printf,1,'echo ".r /astro/users/christensen/code/IDL/HIcubes/analyzeHIcubes.pro" >> analyzeHIcubes.idl.batch'
;printf,1,'echo "analyzeHIcubes,dir,file" >> analyzeHIcubes.idl.batch'
IF keyword_set(dwarfv) THEN BEGIN
    printf,1,'echo "analyzeHIcubes,dir,file,angle = 90,/dwarf,/physical_coord" >> analyzeHIcubes.idl.batch'
    printf,1,'echo "analyzeHIcubes,dir,file,angle = 45,/dwarf,/physical_coord" >> analyzeHIcubes.idl.batch'
ENDIF ELSE BEGIN
    printf,1,'echo "analyzeHIcubes,dir,file,angle = 90,/physical_coord" >> analyzeHIcubes.idl.batch'
    printf,1,'echo "analyzeHIcubes,dir,file,angle = 45,/physical_coord" >> analyzeHIcubes.idl.batch'
ENDELSE
printf,1,'idl analyzeHIcubes.idl.batch'
close,1
spawn,'chmod a+x con_script.analyzeHIcubes'
end
