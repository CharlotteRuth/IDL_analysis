pro testLTE,yH,Radr_HII,TotrForm_H2,Coll_HI,Coll_e_H2,Coll_H_H2

  Phot_HI = 1e-100
  Phot_H2 = 1e-100  

  fHII = (Radr_HII)/(Phot_HI + Coll_HI) ;//rcirrH2 + rpirrH2 * rye;
  fHI = 2.0*(TotrForm_H2)/(Coll_e_H2 + Coll_H_H2 + Phot_H2) ;//rcirrHI + rpirrHI * rye;
  rfH = 1 / ( 1 + fHII * (1 + fHI) ) ; 
  yHII = yH * rfH               ;
  yHI = yH* rfH * fHII          ;
  yH2 = (yH - yHI - yHII)/2.0
  print,fHII,fHI,rfH
  print,yHI,yH2,yHII,yHI+2.0*yH2+yHII,yH

end
