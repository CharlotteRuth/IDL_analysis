FUNCTION TIPSYUNITS,pfile,scale = scale,SILENT = silent,VERBOSE = verbose,ICfile = ICfile
;takes an input gasoline .param file and returns the time unit,
;length of one time step and length of simulation.

;'units=[-1.,-1.,-1.,-1.,-1.,-1,-1.,-1.] ;for returning
;units=[massunit, lengthunit, timeunit,timestep,timetotal, vunit,rhounit, h]
units=CREATE_STRUCT('massunit',-1.0, 'lengthunit',-1.0, 'timeunit',-1.0,'timestep',-1.0,'timetotal',-1.0, 'vunit',-1.0, 'rhounit',-1.0, 'h',-1.0,'gasparmass',-1.0, 'darkparmass',-1.0,'istarmass',-1.0)
cm_per_kpc = 3.08568021d21
gm_per_msol = 1.98892d33
amu_per_gm = 6.022142d23
gm_per_H = 1.673534d-24
H_per_gm = 5.9753790e+23
IF NOT keyword_set(a) THEN scale = 1
READCOL, pfile,param,equals,value,FORMAT='A,A,A',/SILENT

;only look through items that are not commented and are properly formatted
a=STRMATCH(equals,'=')
b=STRMATCH(param,'*#*')
ind=WHERE(a EQ 1 AND b EQ 0)
param=param[ind]
value=value[ind]

massind=WHERE(param EQ 'dMsolUnit',nmassind)
IF (nmassind EQ 1) THEN BEGIN
    massunit=DOUBLE(value[massind[0]])
ENDIF ELSE BEGIN
    print, "No dMsolUnit"
    return,massind
ENDELSE

lengthind=WHERE(param EQ 'dKpcUnit',nlengthind)
IF (nlengthind EQ 1) THEN BEGIN
    lengthunit=DOUBLE(value[lengthind[0]])
ENDIF ELSE BEGIN
    print, "No dKpcUnit"
    return,lengthind
ENDELSE

deltaind=WHERE(param EQ 'dDelta',ndeltaind)
IF (ndeltaind EQ 1) THEN BEGIN
    ddelta=DOUBLE(value[deltaind[0]])
ENDIF ELSE BEGIN
    print, "No dDelta"
    ddelta = 1.0
;    return,deltaind
ENDELSE

stepind=WHERE(param EQ 'nSteps',nstepind)
IF (nstepind EQ 1) THEN BEGIN
    nstep=DOUBLE(value[stepind[0]])
ENDIF ELSE BEGIN
    print, "No nSteps"
    return,stepind
ENDELSE

lambdaind=WHERE(param EQ 'dLambda',lambda)
IF (lambda EQ 1) THEN BEGIN
    lambda=DOUBLE(value[lambdaind[0]])
ENDIF ELSE BEGIN
    print, "No Lambda"
    lambda = 0
;    return,lambdaind
ENDELSE

omegaind=WHERE(param EQ 'dOmega0',omega)
IF (omega EQ 1) THEN BEGIN
    omega=DOUBLE(value[omegaind[0]])
ENDIF ELSE BEGIN
    print, "No Omega"
    omega = 0
;    return,omegaind
ENDELSE

hubbleind=WHERE(param EQ 'dHubble0',hubble)
IF (hubble EQ 1) THEN BEGIN
    hubble=DOUBLE(value[hubbleind[0]])
ENDIF ELSE BEGIN
    print, "No Hubble0"
    hubble = 0
;    return,hubbleind
ENDELSE

imassind=WHERE(param EQ 'dInitStarMass',istarmass)
IF (istarmass EQ 1) THEN BEGIN
    istarmass=DOUBLE(value[imassind[0]])*massunit
    gasparmass = istarmass/0.4
ENDIF ELSE BEGIN
    print, "No InitStarMass"
    istarmass = 0
    gasparmass = 0
ENDELSE

inFileind=WHERE(param EQ 'achInFile',inFile)
IF (inFileind GE 0 AND KEYWORD_SET (ICFILE)) THEN BEGIN
    inFile=STRING(value[inFileind[0]])
    lasti = STRPOS(pfile,'/',/REVERSE_SEARCH)
    path = STRMID(pfile,0,lasti + 1)
    rtipsy,path+inFile,head,g,d,s
    if istarmass eq 0 AND keyword_set(s) then istarmass = MAX(s.mass)*massunit
    if istarmass eq 0 then gasparmass = MAX(g.mass)*massunit   
    darkparmass = MAX(d.mass)*massunit
ENDIF ELSE BEGIN
    print, "No IC"
    darkparmass = 0
ENDELSE

rhounit=massunit * gm_per_msol * H_per_gm/lengthunit^3/cm_per_kpc^3
IF (FLOAT(lambda) eq 0.7) then h = 0.7 ELSE BEGIN
    IF (FLOAT(lambda) eq 0.76) then h = 0.73 ELSE h = 0.71
ENDELSE
vunit = 100.0* h * (lengthunit / 1000.0)/2.894405
timeunit=SQRT((lengthunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
timestep=timeunit*ddelta
timetotal=timestep*nstep
units.massunit = massunit
units.lengthunit = lengthunit*scale
units.timeunit = timeunit
units.timestep = timestep
units.timetotal = timetotal
units.vunit = vunit*scale
units.rhounit = rhounit/scale^3
units.h = h
units.gasparmass = gasparmass
units.darkparmass = darkparmass
units.istarmass = istarmass

;units=[massunit, lengthunit, timeunit,timestep,timetotal, vunit, rhounit, h, gasparmass, darkparmass, istarmass]
IF (KEYWORD_SET(verbose)) THEN BEGIN
    print, "Mass Unit [Msol]: ",massunit
    print, "Length Unit [kpc]: ",lengthunit
    print, "Time Unit [yrs]: ", timeunit
    print, "Time Step [yrs]: ", timestep
    print, "Sim. Run Time [yrs]: ", timetotal
    print, "Velocity Unit [km s^-1]: ",vunit
    print, "Density Unit [atoms cm^-3]: ",rhounit
    print, "h0 [km s^-1 Mpc^-1]: ",h
    IF (inFileind GE 0 AND KEYWORD_SET (ICFILE)) THEN print, "Init. number of dark particles: ",head.ndark
    IF (inFileind GE 0 AND KEYWORD_SET (ICFILE)) THEN print, "Init. number of gas particles: ",head.ngas
    IF (inFileind GE 0 AND KEYWORD_SET (ICFILE)) THEN print, "Dark particle mass [Msol]: ",MINMAX(d.mass)*massunit
    IF (inFileind GE 0 AND KEYWORD_SET (ICFILE)) THEN print, "Gas particle mass [Msol]: ",MINMAX(g.mass)*massunit
    print, "Stellar mass (init) [Msol]: ",istarmass
ENDIF
return, units
END
