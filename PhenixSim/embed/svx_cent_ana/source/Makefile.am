
AUTOMAKE_OPTIONS = foreign

BUILT_SOURCES = \
  testexternals.C


INCLUDES = -I. -I$(srcdir) -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = libsvxcentana.la

noinst_HEADERS = svxcentanaLinkDef.h 

include_HEADERS = svxana.h embedana.h embeddisp.h SvxCentralTrackRecoDebug.h svxcentclsana.h svxstripdisp.h svxcntclusdisp.h svxcnttrkdisp.h SvxCheckCompactCNT.h ReadNtpCntTrk.h AnaNtpTrkPair.h AnaCompactCNT.h SvxanaPair.h svxanamini.h DebugTest.h SvxRpQA.h SvxanaPairDisp.h SvxCheckCluster.h

libsvxcentana_la_SOURCES = svxana.cc embedana.cc embeddisp.cc SvxCentralTrackRecoDebug.C svxcentclsana.cc svxstripdisp.cc CentAnaDict.C svxcntclusdisp.cc svxcnttrkdisp.cc SvxCheckCompactCNT.C ReadNtpCntTrk.C AnaNtpTrkPair.C AnaCompactCNT.C SvxanaPair.C svxanamini.cc DebugTest.C SvxRpQA.C SvxanaPairDisp.C SvxCheckCluster.C

libsvxcentana_la_LIBADD = -L$(OFFLINE_MAIN)/lib -lphool -lfun4all -luti -lCNT -lsvx -lreactionplane -lsimreco -lcompactCNT -lrecal


noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = \
  testexternals.C

testexternals_LDADD = \
  libsvxcentana.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@


CentAnaDict.C: svxana.h embedana.h embeddisp.h SvxCentralTrackRecoDebug.h svxcentclsana.h svxstripdisp.h svxcntclusdisp.h svxcnttrkdisp.h SvxCheckCompactCNT.h ReadNtpCntTrk.h AnaNtpTrkPair.h AnaCompactCNT.h  SvxanaPair.h svxanamini.h DebugTest.h SvxRpQA.h SvxanaPairDisp.h SvxCheckCluster.h svxcentanaLinkDef.h
	rootcint -f $@ -c $(INCLUDES) $^

# Rule for generating CINT dictionaries from class headers.
%Dict.C %Dict.h: %.hh
	rootcint -f *Dict* -c $(INCLUDES) $*.hh

