FUNCTION TIMEUNITS,pfile,SILENT=silent
;takes an input gasoline .param file and returns the time unit,
;length of one time step and length of simulation.

timevals=[-1.,-1.,-1.] ;for returning


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
    return,deltaind
ENDELSE

stepind=WHERE(param EQ 'nSteps',nstepind)
IF (nstepind EQ 1) THEN BEGIN
    nstep=DOUBLE(value[stepind[0]])
ENDIF ELSE BEGIN
    print, "No nSteps"
    return,stepind
ENDELSE

timeunit=SQRT((lengthunit*3.086d21)^3/(6.67d-8*massunit*1.99d33))/(3600.*24.*365.24)
timestep=timeunit*ddelta
timetotal=timestep*nstep
timevals=[timeunit,timestep,timetotal]
IF (KEYWORD_SET(silent) EQ 0) THEN BEGIN
    print, "Time Unit [yrs]: ", timeunit
    print, "Time Step [yrs]: ", timestep
    print, "Sim. Run Time [yrs]: ", timetotal
ENDIF
return, timevals
END
