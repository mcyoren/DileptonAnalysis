#include "Run14AuAuLeptonCombyReco.h"

bool Run14AuAuLeptonCombyReco::IsCentralSupportCut(const float theta0, const float bbcVertex)
{
    
    if(
        !(  fabs(theta0) < 100 &&
            (
                (
                    bbcVertex > 0 &&
                    (
                        bbcVertex - 250.*tan(theta0 - 3.1416/2.) > 2 ||
                        bbcVertex - 200.*tan(theta0 - 3.1416/2.) < -2
                    )
                ) ||
                (
                    bbcVertex < 0 &&
                    (
                        bbcVertex - 250.*tan(theta0 - 3.1416/2.) < -2 ||
                        bbcVertex - 200.*tan(theta0 - 3.1416/2.) > 2
                    )
                )
            )
         )
    )
        return true;

    return false;
}
