# by default, compile only
default: bin/mesh_slicer

clean:
	$(RM) -f bin/mesh_slicer

# compile the slicing executable
CGAL_INSTALL_PREFIX ?= /opt/homebrew
CGAL_INC_DIR = $(CGAL_INSTALL_PREFIX)/include
CGAL_LIB_DIR = $(CGAL_INSTALL_PREFIX)/lib
CXXFLAGS = -O2 -std=c++17
LDFLAGS = -Wl,-rpath,$(CGAL_LIB_DIR)
LDLIBS = -lgmp -lmpfr

bin/mesh_slicer: src/mesh_slicer.cpp
	$(CXX) $(CXXFLAGS) -I$(CGAL_INC_DIR) $< -o $@ -L$(CGAL_LIB_DIR) $(LDFLAGS) $(LDLIBS)

# convert STL to OFF using PyMeshLab
PYTHON ?= python3

%.off: %.stl
	$(PYTHON) bin/removeDuplicatedVertex.py $< $@

# ======= FORMAT FOR TARGETS =======
#output: executable meshfile centrelinefile
#	./$^ $@   # this works due to the order of arguments and above dependencies

# demo with cylinder
.PHONY: demo
demo: demo/slices.dat
demo/slices.dat: bin/mesh_slicer demo/cylinder.off demo/centreline.dat
	./$^ $@
