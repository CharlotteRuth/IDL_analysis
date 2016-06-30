;
;
;
;  use peter's disk fitting routine to plot scaleheights of old and
;  young stars
;
;
;
;


pro scaleheight, s

bin = 0.2
ls -
r = sqrt(s.x^2+s.y^2)

ind = where(r gt 5)

im = better_hist2d(s[ind].x,s[ind].z,s[ind].mass, binsize1=bin, binsize2=bin, $
                   obin1=obin1, obin2=obin2)

zprof = total(im,1)

stop

end
