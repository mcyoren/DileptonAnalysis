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
      ZVTX,
      PHI1,
      PHI2,
      PHI3,
      THE1,
      THE2,
      THE3,
      WEIGHT,
      LAST_DOUBLE  // special...always use this
    };

  enum Run14AuAuLeptonCombyEnumIntegerFields
    {
      PTYPE,
      HADRON_REJECT,
      HIT_ASSOC,
      CONV_REJECT,
      GHOST,
      SECTOR,
      YSECT,
      ZSECT,
      CENTRALITY,
      BBCQ,
      LAST_INTEGER  // special...always use this
    };
};

#endif
