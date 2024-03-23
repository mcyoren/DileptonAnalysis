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
      //DCA,
      //SDCA,
      PHI0,
      DCAX,
      DCAY,
      KEFF,
      LAST_DOUBLE  // special...always use this
    };

  enum Run14AuAuLeptonCombyEnumIntegerFields
    {
      MATCH,
      GHOST,
      PTYPE,
      LAST_INTEGER  // special...always use this
    };
};

#endif
