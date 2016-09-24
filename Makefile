PLUGIN_NAME = livshitz_rudy_2009

RTXI_INCLUDES=/usr/local/lib/rtxi_includes

HEADERS = LivR2009_Model.h

SOURCES = LivR2009_Model.cpp \
          ${RTXI_INCLUDES}/powfast.cpp

LIBS = 

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
