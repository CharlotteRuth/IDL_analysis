;This function returns a tipsy particle whoes position has been
;shifted over by the center vector
FUNCTION center_particle,particle,center

particle.x = particle.x - center[0]
particle.y = particle.y - center[1]
particle.z = particle.z - center[2]

RETURN,particle
END
