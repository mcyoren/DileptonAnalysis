AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  field_map.root

lib_LTLIBRARIES = \
  libRun14AuAuLeptonEvent.la \
  libRun14AuAuLeptonComby.la 

noinst_HEADERS = \
  MyEventLinkDef.h \
  MyEvent.h \
  MyEventConstants.h \
  Run14AuAuLeptonCombyLinkDef.h

include_HEADERS = \
  Run14AuAuLeptonCombyEnum.h \
  Run14AuAuLeptonCombyConstants.h \
  Run14AuAuLeptonCombyReco.h \
  Declarations.h \
  Run14AuAuLeptonCombyHistos.h \
  Run14AuAuLeptonCombyCutter.h 

libRun14AuAuLeptonEvent_la_SOURCES = \
  MyEvent.C \
  MyEvent_Dict.C

libRun14AuAuLeptonComby_la_SOURCES = \
  Run14AuAuLeptonCombyReco.C \
  IsCentralSupportCut.C \
  QACuts.C \
  Run14AuAuLeptonCombyHistos.C \
  Run14AuAuLeptonCombyCutter.C \
  Run14AuAuLeptonComby_Dict.C

libRun14AuAuLeptonEvent_la_LIBADD = \
$(ROOT_LIBS)  

libRun14AuAuLeptonComby_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lTOAD \
  -lfun4allfuncs \
  -lMWGpico \
  -lCNT \
  -lUltraLight \
  -lCabanaBoy

# Rule for generating CINT dictionaries from class headers.
MyEvent_Dict.C: \
  MyEvent.h \
  MyEventConstants.h \
  MyEventLinkDef.h 
	  rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

# Rule for generating CINT dictionaries from class headers.
Run14AuAuLeptonComby_Dict.C: \
  Run14AuAuLeptonCombyReco.h \
  Declarations.h \
  Run14AuAuLeptonCombyEnum.h \
  Run14AuAuLeptonCombyConstants.h \
  Run14AuAuLeptonCombyCutter.h \
  Run14AuAuLeptonCombyHistos.h \
  Run14AuAuLeptonCombyLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

