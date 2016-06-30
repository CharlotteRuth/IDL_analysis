;This function reads in skid output .stat files and puts them in an array.

FUNCTION read_stat, filename

;test the file to make sure there are halos
test=FILE_TEST(filename,/ZERO)

;if file has non-zero length
IF (test EQ 0) THEN BEGIN

    readcol, filename, group, members, mtot, mgas, mstar, vcmax, vchalf, vcouter, r_vcmax, r_mhalf, r_outer, vdisp, x, y, z, vx, vy, vz, xbound, ybound, zbound,FORMAT='(L,L,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F)',/SILENT

    ngroups=N_ELEMENTS(group)

ENDIF ELSE BEGIN
    ;define the variables so halo still contains real data.
    group=0l & members=0l
    mtot=0. & mgas=0. & mstar=0. & vcmax=0. & vchalf=0. & vcouter=0.
    r_vcmax=0. & r_mhalf=0. & r_outer=0. & vdisp=0. & x=0. & y=0. 
    z=0. & vx=0. & vy=0. & vz=0. & xbound=0. & ybound=0. & zbound=0.
ENDELSE
    
;define halo structure to hold all the data.
halo={group:group, members:members, mtot:mtot, mgas:mgas, mstar:mstar, vcmax:vcmax, vchalf:vchalf, vcouter:vcouter, r_vcmax:r_vcmax, r_mhalf:r_mhalf, r_outer:r_outer, vdisp:vdisp, x:x, y:y, z:z, vx:vx, vy:vy, vz:vz, xbound:xbound, ybound:ybound, zbound:zbound}

RETURN, halo

END
