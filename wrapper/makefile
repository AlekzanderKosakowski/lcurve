# ------------------------------------------------------------
# Detect the Python interpreter and its include directory
# ------------------------------------------------------------
PYTHON      := python3
PY_VER      := $(shell $(PYTHON) -c "import sys; \
                      print('{}.{}'.format(sys.version_info[0], sys.version_info[1]))")
PY_INC      := -I$(shell $(PYTHON) -c "import sysconfig; \
                      print(sysconfig.get_paths()['include'])")

# ------------------------------------------------------------
# Legacy project include directories (renamed)
# ------------------------------------------------------------
CPP_INCS    := \
    -I../cpp-binary/include \
    -I../cpp-roche/include \
    -I../cpp-lcurve/include \
    -I../cpp-colly/include \
    -I../cpp-subs/include

# ------------------------------------------------------------
# Local pybind11 headers (cloned into external/)
# ------------------------------------------------------------
PYB_INCS    := -I../external/pybind11/include

# ------------------------------------------------------------
# PGPLOT header – absolute system location
# ------------------------------------------------------------
PGPLOT_INCS := -I/usr/local/src/pgplot/include

# ------------------------------------------------------------
# Combine all include flags (each entry already has the leading -I)
# ------------------------------------------------------------
INCLUDE_FLAGS = $(PY_INC) $(CPP_INCS) $(PYB_INCS) $(PGPLOT_INCS)

# ------------------------------------------------------------
# Compiler and flags
# ------------------------------------------------------------
CXX         := g++
CXXFLAGS    := -O3 -Wall -std=c++14 -fPIC $(INCLUDE_FLAGS)


# Add the PCRE‑C++ libraries (order matters: pcrecpp depends on pcre)
LDLIBS     := -L/usr/local/lib -lcsla -lpcrecpp -lpcre -lgomp

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------
TARGET_SO   := lcurve_wrapper.cpython-$(PY_VER)-x86_64-linux-gnu.so
STATIC_LIB  := ../lib/liblegacy.a
SRC         := wrapper.cpp

# ------------------------------------------------------------
# Build rules
# ------------------------------------------------------------
all: $(TARGET_SO)

$(TARGET_SO): $(SRC) $(STATIC_LIB)
	@echo "=== COMPILING $(TARGET_SO) ==="
	@echo "$(CXX) $(CXXFLAGS) -shared $(SRC) $(STATIC_LIB) -o $(TARGET_SO)"
	$(CXX) $(CXXFLAGS) -shared $(SRC) $(STATIC_LIB) $(LDLIBS) -o $@

# Optional: rebuild the static archive if any .o files changed
$(STATIC_LIB):
	@echo "=== REBUILDING static archive liblegacy.a ==="
	ar rcs $(STATIC_LIB) \
	    ../cpp-binary/src/*.o \
	    ../cpp-roche/src/*.o \
	    ../cpp-lcurve/src/*.o \
	    ../cpp-colly/src/*.o \
	    ../cpp-subs/src/*.o

clean:
	rm -f $(TARGET_SO)

.PHONY: all clean
