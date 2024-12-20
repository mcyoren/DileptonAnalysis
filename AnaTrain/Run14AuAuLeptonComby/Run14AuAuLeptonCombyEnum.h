#ifndef __RUN14AUAULEPTONCOMBYENUM_H__
#define __RUN14AUAULEPTONCOMBYENUM_H__

class Run14AuAuLeptonCombyEnum {

public:

    enum Run14AuAuLeptonCombyEnumDoubleFields {
      ALPHA,
      PHI,
      ZED,
      CRKPHI,
      CRKZED,
      DCAX,
      DCAY,
      PSI,
      ZVTX,
      PHI1,
      PHI2,
      PHI3,
      THE1,
      THE2,
      THE3,
      LAST_DOUBLE  // special...always use this
    };

  enum Run14AuAuLeptonCombyEnumIntegerFields
    {
      MATCH,
      GHOST,
      PTYPE,
      ID1,
      ID2,
      ID3,
      CENTR,
      LAST_INTEGER  // special...always use this
    };
};

#endif
