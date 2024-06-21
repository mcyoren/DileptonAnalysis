AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

mydatadir = $(datadir)/$(PACKAGE_INCLUDE_DIR)
dist_mydata_DATA = \
  field_map.root

lib_LTLIBRARIES = \
  libRun14AuAuLeptonEvent.la \
  libsvxcentana.la

noinst_HEADERS = \
  MyEventLinkDef.h \
  MyEvent.h \
  MyEventConstants.h 

include_HEADERS = \
  embedana.h 

libRun14AuAuLeptonEvent_la_SOURCES = \
  MyEvent.C \
  MyEvent_Dict.C

libsvxcentana_la_SOURCES = \
  embedana.cc \
  embedana_Dict.C

libRun14AuAuLeptonEvent_la_LIBADD = \
$(ROOT_LIBS)  

libsvxcentana_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lTOAD \
  -lphool \
  -lfun4all \
  -luti \
  -lCNT \
  -lsvx \
  -lreactionplane \
  -lsimreco \
  -lcompactCNT \
  -lrecal

# Rule for generating CINT dictionaries from class headers.
MyEvent_Dict.C: \
  MyEvent.h \
  MyEventConstants.h \
  MyEventLinkDef.h 
	  rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

# Rule for generating CINT dictionaries from class headers.                                                                                                                                     
embedana_Dict.C: \
  embedana.h \
  embedanaLinkDef.h 
	  rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^


clean-local:
	rm -f *Dict*

