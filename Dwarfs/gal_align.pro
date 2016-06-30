function gal_align,h,g,d,s, RBOX = rbox,  LUNIT = lunit, verbose=verbose

; NAME:
;   bw_rotgal
;
; PURPOSE:
;     To create a TIPSY file of some or all of a simulation snapshot,
;     rotated such that the angular momentum vector of the primary
;     galaxy is oriented along the z-direction, with an optional
;     angular offset.  
;
; CALLING SEQUENCE:
;  IDL> bw_rotgal, infile = infile
;
; REQUIRED INPUTS:
;     infile = Name of TIPSY snapshot file
;      e.g. 'MW1.1024g1bwK.00512'
;
; OPTIONAL INPUTS:
;
;     rbox = In comoving kpc, the halfwidth of a cube within which particles are
;     included in the angular momentum calculation that will orient
;     the disk.  Generally want
;     this to be small.  Bars can screw you up, so always visually
;     check output with tipsy. Default = 1 comoving kpc.
;
;     outbox = In physical units, the DIAMETER of the box you would
;     like to output in the rotated coordinate system. Default =
;     entire box. If you want the whole box, do not set outbox.
;                                                                
;
;     angle = Sets an angular offset in the rotated coordinate
;     system. Angle of 0 would yield an output box with the galaxy
;     angular momentum vector in the z-direction (edge-on), and an
;     angle of 90 would yield a face on system.  Default - 0.  Only
;                                                          accepts
;                                                          angle > 0.
;
;     lunit = Comoving length unit of the simulation.  Used to convert the rbox
;     and inbox to simulation coordinates. Default = 2.85714e4 kpc 
;
;     type = type of particle used to calculate angular momentum of galaxy.
;     Acceptable inputs are: gas, star, or baryon. Default = gas.
;
;     haloid = The AMIGA grp number assigned to the halo that is being 
;     rotated.  Default = 1 (typically the most massive galaxy in the 
;     simulation. 
;
;     multiplefile = if you wish to orient and cut multiple halos from 
;     one simulation, this file contains the optional inputs for each 
;     halo.  This allows rtipsy to be called just once for the simulation, 
;     drastically reducing runtime.  In this case, the input file MUST 
;     contain a column for EACH of the optional inputs.
;
; OUTPUTS:
;     A rotated (sub)set of the input snapshot in TIPSY format in
;     simulation units.  
;       infile+".size.angle.haloid.rot", where "size" is outbox and "angle"
;       is the number of degrees between disk angular momentum vector
;       and z vector, and haloid is the grp number of the halo.  
;     e.g. MW1.1024g1bwK.00512.025.00.01.rot is a 25 kpc diamter box 
;     around the MW main halo with the angular momentum vector 
;     aligned with the z-direction.
;
;     The size 999 corresponds to a rotation of the whole box.
;     e.g. MW1.1024g1bwK.00512.999.00.rot
;
;     A markfile that contains the indices of all particles in the
;     output box.  The markfile has the format:
;       infile+".size.angle.rot.mark"
;     A markfile is not output for a whole box rotation
;
; EXTERNAL ROUTINES:
;     rtipsy
;     wtipsy
;     az_norm
;     az_dotp
;     readcol
;     
; EXTERNAL DATA:

; COMMENTS:
;   LIMITATIONS = Only accepts angle > 0.  CURRENTLY DOESN"T ACCEPT
;   NON-ZERO ANGLE.  DOES NOT YET OUTPUT A MARKFILE.  The "angle" 
;   input is thus currently ignored. 
;   
;   Currently the code relies on the halo finder being amiga, the halo
;   finder naming convention being infile+.amiga.gtp, and the first
;   element of the .gtp file being the primary halo of the sim box
;
;   It can be helpful to output only a subset of the snapshot in
;   rotation form for purposes of image creation and running the
;   galaxy scope.
;
; REVISION HISTORY:
;   2006-October-28 Alyson Brooks (added ability to choose the haloid 
; 		number, and made "outbox" work) 
;   2006-October-13 Beth Willman, version 2.0
;   2006-August     Adi Zolotov, version 1.0 writes rotation matric


IF n_elements(lunit) eq 0 then lunit=2.857d4 ; kpc per sim unit
; Put rbox into simulation units
rbox = rbox/lunit

;Select box of particles to use
;dist_x=g.X-halo[haloind].x
;dist_y=g.Y-halo[haloind].y
;dist_z=g.Z-halo[haloind].z
dist_x=g.X
dist_y=g.Y
dist_z=g.Z
pcuts=where((abs(dist_x) lt rbox) and $
            (abs(dist_y) lt rbox) and $
            (abs(dist_z) lt rbox),ncuts)

print, "Number of gas particles in L box:", ncuts
IF (N_ELEMENTS(pcuts) GT 5) then begin 
    box=g[pcuts]
    print,'Enough gas particles to continue'
;    x=box.x-halo[haloind].x
;    y=box.y-halo[haloind].y
;    z=box.z-halo[haloind].z
    x=box.x
    y=box.y
    z=box.z
    px=(box.vx)*box.mass
    py=(box.vx)*box.mass
    pz=(box.vz)*box.mass
    px=(box.vx)*box.mass
    py=(box.vy)*box.mass
    pz=(box.vz)*box.mass

    L_box= dblarr(3)
    for ii=0L,n_elements(box)-1L do L_box $
       = L_box +crossp([x[ii],y[ii],z[ii]],[px[ii],py[ii],pz[ii]])
    
    zhat=L_box/az_norm(L_box)

;begplot, 'rotation.ps',xsize=7.0,ysize=7.0
;xrange=[-1,1]*rbox
;yrange=[-1,1]*rbox
;window,1
;plot,x,y,psym=3,title='before',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[0]],[0.0,zhat[1]]
;stop
;plot,x,z,psym=3,title='before',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[0]],[0.0,zhat[2]]
;stop
;plot,y,z,psym=3,title='before',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[1]],[0.0,zhat[2]]
;stop
;plot,g.x,g.y,psym=3,title='before'

    r_nth=[x[1],y[1],z[1]]
    xtilda=r_nth-zhat*az_dotp(zhat,r_nth)
    xhat=xtilda/az_norm(xtilda)
    yhat=crossp(zhat,xhat)

    IF KEYWORD_SET(verbose) then begin
;print out to make sure everything is ok
;        splog,xhat
;        splog,yhat
;        splog,zhat

 ;       splog, az_dotp(xhat,xhat)
 ;       splog, az_dotp(xhat,yhat)
 ;       splog, az_dotp(xhat,zhat)
 ;       splog, crossp(xhat,yhat)-zhat
 ;       splog, crossp(yhat,zhat)-xhat
 ;       splog, crossp(zhat,xhat)-yhat
        print,xhat
        print,yhat
        print,zhat

        print, az_dotp(xhat,xhat)
        print, az_dotp(xhat,yhat)
        print, az_dotp(xhat,zhat)
        print, crossp(xhat,yhat)-zhat
        print, crossp(yhat,zhat)-xhat
        print, crossp(zhat,xhat)-yhat
   ENDIF

        axis=dblarr(3,3)
        axis[*,0]=xhat
        axis[*,1]=yhat
        axis[*,2]=zhat

;rotate positions of gas

        gx = (axis##TEMPORARY([[g.x],[g.y],[g.z]]))[*,0]
        gy = (axis##TEMPORARY([[g.x],[g.y],[g.z]]))[*,1]
        g.z = (axis##TEMPORARY([[g.x],[g.y],[g.z]]))[*,2]
        g.x = gx
        g.y = gy
  
        gvx = (axis##TEMPORARY([[g.vx],[g.vy],[g.vz]]))[*,0]
        gvy = (axis##TEMPORARY([[g.vx],[g.vy],[g.vz]]))[*,1]
        g.vz = (axis##TEMPORARY([[g.vx],[g.vy],[g.vz]]))[*,2]
        g.vx = gvx
        g.vy = gvy

        dx = (axis##TEMPORARY([[d.x],[d.y],[d.z]]))[*,0]
        dy = (axis##TEMPORARY([[d.x],[d.y],[d.z]]))[*,1]
        d.z = (axis##TEMPORARY([[d.x],[d.y],[d.z]]))[*,2]
        d.x = dx
        d.y = dy
   
        dvx = (axis##TEMPORARY([[d.vx],[d.vy],[d.vz]]))[*,0]
        dvy = (axis##TEMPORARY([[d.vx],[d.vy],[d.vz]]))[*,1]
        d.vz = (axis##TEMPORARY([[d.vx],[d.vy],[d.vz]]))[*,2]
        d.vx = dvx
        d.vy = dvy

        sx = (axis##TEMPORARY([[s.x],[s.y],[s.z]]))[*,0]
        sy = (axis##TEMPORARY([[s.x],[s.y],[s.z]]))[*,1]
        s.z = (axis##TEMPORARY([[s.x],[s.y],[s.z]]))[*,2]
        s.x = sx
        s.y = sy       

        svx = (axis##TEMPORARY([[s.vx],[s.vy],[s.vz]]))[*,0]
        svy = (axis##TEMPORARY([[s.vx],[s.vy],[s.vz]]))[*,1]
        s.vz = (axis##TEMPORARY([[s.vx],[s.vy],[s.vz]]))[*,2]
        s.vx = svx
        s.vy = svy

;window,2
;xrange=[-1,1]*rbox
;yrange=[-1,1]*rbox
;plot,g[pcuts].x,g[pcuts].y,psym=3,title='After',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[0]],[0.0,zhat[1]]
;stop
;plot,g[pcuts].x,g[pcuts].z,psym=3,title='After',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[0]],[0.0,zhat[2]]
;stop
;plot,g[pcuts].y,g[pcuts].z,psym=3,title='After',xrange=xrange,yrange=xrange
;oplot, [0.0,zhat[1]],[0.0,zhat[2]]
;stop
;plot,g.x,g.y,psym = 3
ENDIF
particles = {h:h,g:g,d:d,s:s}
return,particles
END
