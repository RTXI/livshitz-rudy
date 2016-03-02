PLUGIN_NAME = livshitz_rudy_2009

HEADERS = LivR2009_Model.h

SOURCES = LivR2009_Model.cpp \
          include/PowFast.cpp

LIBS = 

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
