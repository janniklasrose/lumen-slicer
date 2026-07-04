# by default, compile only
.PHONY: default
default: bin/mesh_slicer bin/tessellate_slices bin/surface_reconstructor

.PHONY: clean
clean:
	$(RM) bin/mesh_slicer bin/tessellate_slices bin/surface_reconstructor

# compile the executables
CGAL_INSTALL_PREFIX ?= /opt/homebrew
CGAL_INC_DIR = $(CGAL_INSTALL_PREFIX)/include
CGAL_LIB_DIR = $(CGAL_INSTALL_PREFIX)/lib
CXXFLAGS = -O2 -std=c++17
LDFLAGS = -Wl,-rpath,$(CGAL_LIB_DIR)
LDLIBS = -lgmp -lmpfr

bin/%: src/%.cpp
	$(CXX) $(CXXFLAGS) -I$(CGAL_INC_DIR) $< -o $@ -L$(CGAL_LIB_DIR) $(LDFLAGS) $(LDLIBS)

# convert STL to OFF using PyMeshLab
PYTHON ?= python3

%.off: %.stl
	$(PYTHON) bin/removeDuplicatedVertex.py $< $@

.PHONY: demo
demo: demo/slices.dat demo/slices.off demo/reconstructed.off

# ======= FORMAT FOR TARGETS =======
#output: executable meshfile centrelinefile
#	./$^ $@   # this works due to the order of arguments and above dependencies

# convert cylinder to slices
demo/slices.dat: bin/mesh_slicer demo/cylinder.off demo/centreline.dat
	./$^ $@
# tessellate the slices
demo/slices.off: bin/tessellate_slices demo/slices.dat
	./$^ $@
# reconstruct the outer surface from the slices
demo/reconstructed.off: bin/surface_reconstructor demo/slices.dat
	./$^ $@
