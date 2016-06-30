; center the particles around their center of mass

pro cofm, header, g, d, s, POT = pot, cofm_xyz = cofm_xyz, cofmv_xyz = cofmv_xyz, $
          nocenter = nocenter

if(header.ngas eq 0) then begin
    g = replicate({mass: 0.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,$
                   dens:1.,tempg:1.,h : 1. , zmetal : 1., phi : 1.},1)
endif

if(header.ndark eq 0) then begin
    d = replicate({mass: 0.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,$
                   eps: 1.,phi: 1.},1)
endif

if(header.nstar eq 0) then begin
    s = replicate({mass: 0.,x: 1.,y : 1., z:1.,vx:1.,vy:1.,vz:1.,$
                   metals:1.,tform:1.,eps: 1.,phi: 1.},1)
endif


if(keyword_set(pot)) then begin
    g_pot_min_ind = where(g.phi eq min(g.phi))
    g_pot_min = g[g_pot_min_ind].phi

    if(n_elements(g_pot_min_ind) gt 1) then begin
        g_pot_min_ind = g_pot_min_ind[0]
        g_pot_min = g_pot_min[0]
    endif
    
    d_pot_min_ind = where(d.phi eq min(d.phi))
    d_pot_min = d[d_pot_min_ind].phi

    if(n_elements(d_pot_min_ind) gt 1) then begin
        d_pot_min_ind = d_pot_min_ind[0]
        d_pot_min = d_pot_min[0]
    endif

    s_pot_min_ind = where(s.phi eq min(s.phi))
    s_pot_min = s[s_pot_min_ind].phi

    if(n_elements(s_pot_min_ind) gt 1) then begin
        s_pot_min_ind = s_pot_min_ind[0]
        s_pot_min = s_pot_min[0]
    endif

    pot_min = [g_pot_min,d_pot_min,s_pot_min]
    type_min = where(pot_min eq min(pot_min))
    
    type_min = type_min[0]
    
    switch type_min of
        0: begin
            cofmx = g[g_pot_min_ind].x
            cofmy = g[g_pot_min_ind].y
            cofmz = g[g_pot_min_ind].z
        end
        1: begin
            cofmx = d[d_pot_min_ind].x
            cofmy = d[d_pot_min_ind].y
            cofmz = d[d_pot_min_ind].z
        end
        2: begin
            cofmx = s[s_pot_min_ind].x
            cofmy = s[s_pot_min_ind].y
            cofmz = s[s_pot_min_ind].z
        end
    endswitch
    
    if(n_elements(cofmx) ne 1) then begin
        cofmx = cofmx[0]
        cofmy = cofmy[0]
        cofmz = cofmz[0]
    endif


    print, 'potential minimum xyz:', cofmx, cofmy, cofmz

endif else begin

    x_coord = [g.x,d.x,s.x]
    y_coord = [g.y,d.y,s.y]
    z_coord = [g.z,d.z,s.z]

    mass = [g.mass, d.mass, s.mass]
    
    cofmx = total(x_coord*mass)/total(mass)
    cofmy = total(y_coord*mass)/total(mass)
    cofmz = total(z_coord*mass)/total(mass)

    cofmx_gas = total(g.x*g.mass)/total(g.mass)
    cofmy_gas = total(g.y*g.mass)/total(g.mass)
    cofmz_gas = total(g.z*g.mass)/total(g.mass)
    
    cofmx_dark = total(d.x*d.mass)/total(d.mass)
    cofmy_dark = total(d.y*d.mass)/total(d.mass)
    cofmz_dark = total(d.z*d.mass)/total(d.mass)

    cofmx_star = total(s.x*s.mass)/total(s.mass)
    cofmy_star = total(s.y*s.mass)/total(s.mass)
    cofmz_star = total(s.z*s.mass)/total(s.mass)

    print, 'cofm gas xyz:', cofmx_gas, cofmy_gas, cofmz_gas
    print, 'cofm dark xyz:', cofmx_dark, cofmy_dark, cofmz_dark
    print, 'cofm star xyz:', cofmx_star, cofmy_star, cofmz_star
    print, 'cofm total xyz:', cofmx, cofmy, cofmz
endelse
    
; calculate the center of mass velocity

x_vel = [g.vx,d.vx,s.vx]
y_vel = [g.vy,d.vy,s.vy]
z_vel = [g.vz,d.vz,s.vz]

mass = [g.mass, d.mass, s.mass]

cofmvx = total(x_vel*mass)/total(mass)
cofmvy = total(y_vel*mass)/total(mass)
cofmvz = total(z_vel*mass)/total(mass)

print, 'cofm vel xyz: ', cofmvx, cofmvy, cofmvz

if(keyword_set(nocenter) eq 0) then begin

; center the particles and set velocities in cofm frame

    g.x = g.x - cofmx
    g.y = g.y - cofmy
    g.z = g.z - cofmz

    s.x = s.x - cofmx
    s.y = s.y - cofmy
    s.z = s.z - cofmz

    d.x = d.x - cofmx
    d.y = d.y - cofmy
    d.z = d.z - cofmz

    g.vx = g.vx - cofmvx
    g.vy = g.vy - cofmvy
    g.vz = g.vz - cofmvz

    s.vx = s.vx - cofmvx
    s.vy = s.vy - cofmvy
    s.vz = s.vz - cofmvz

    d.vx = d.vx - cofmvx
    d.vy = d.vy - cofmvy
    d.vz = d.vz - cofmvz
endif

if(keyword_set(cofm_xyz)) then begin
    cofm_xyz = fltarr(3)
    cofm_xyz[0] = cofmx
    cofm_xyz[1] = cofmy
    cofm_xyz[2] = cofmz
endif

if(keyword_set(cofmv_xyz)) then begin
    cofmv_xyz = fltarr(3)
    cofmv_xyz[0] = cofmvx
    cofmv_xyz[1] = cofmvy
    cofmv_xyz[2] = cofmvz
endif

end
