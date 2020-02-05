#include "global.h"
#include "TCut.h"

void readCut(TCut &cut1, TCut &cut2, TCut &cut3, TCut &cut4, TCut &cut5, TCut &cut6, TCut &cut7, TCut &cut8, TCut &cut9, TCut &cut10, Int_t iPol = -1);

void readCut(TCut cut[], Int_t iPol = -1) {

  cut[0] = "";
  TCut cut1,cut2,cut3,cut4,cut5,cut6,cut7,cut8,cut9,cut10;
  readCut(cut1,cut2,cut3,cut4,cut5,cut6,cut7,cut8,cut9,cut10,iPol);

  cut[ 1]  = cut1;
  cut[ 2]  = cut2;
  cut[ 3]  = cut3;
  cut[ 4]  = cut4;
  cut[ 5]  = cut5;
  cut[ 6]  = cut6;
  cut[ 7]  = cut7;
  cut[ 8]  = cut8;
  cut[ 9]  = cut9;
  cut[10]  = cut10;

}


void readCut(TCut &cut1, TCut &cut2, TCut &cut3, TCut &cut4, TCut &cut5, TCut &cut6, TCut &cut7, TCut &cut8, TCut &cut9, TCut &cut10, Int_t iPol) {

  if (iPol == -1) {
#ifndef __electron__
    cut1  = "leptype==13&&abs(mz-91.2)<10" ;
    cut2  = "npfosc1>3&&npfosc2>3" ;
    cut3  = "evis+elep1+elep2>300" ;      
    cut4  = "bmax1>0.66" ;
    cut5  = "abs(cosz)<0.9" ;
    cut6  = "mhnew>110&&mhnew<150" ;    
    cut7  = "mhnew>110" ;
    cut8  = "mhnew>110" ;
    cut9  = "mhnew>110" ;
    cut10 = "mhnew>110" ;
    //    cut10 = "nhbb==1" ;        
#else
    cut1  = "leptype==11&&abs(mz-91.2)<10" ;
    cut2  = "npfosc1>3&&npfosc2>3" ;
    cut3  = "evis+elep1+elep2>300" ;      
    cut4  = "bmax1>0.66" ;
    cut5  = "abs(cosz)<0.9" ;
    cut6  = "mhnew>110&&mhnew<150" ;    
    cut7  = "mhnew>110" ;
    cut8  = "mhnew>110" ;
    cut9  = "mhnew>110" ;
    cut10 = "mhnew>110" ;    
    //    cut10 = "nhbb==1" ;        
#endif
  }
  else {
  }


}
