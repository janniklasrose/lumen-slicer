# by default, compile only
default: bin/mesh_slicer

clean:
	$(RM) -f bin/mesh_slicer

# compile the slicing executable
CXXFLAGS = -O2 -std=c++11
bin/mesh_slicer: src/mesh_slicer.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ -lCGAL -lgmp

# convert STL to OFF using MeshLab
%.off: %.stl
	meshlabserver -i $< -o $@ -s bin/removeDuplicatedVertex.mlx

# ======= FORMAT FOR TARGETS =======
#output: executable meshfile centrelinefile
#	./$^ $@   # this works due to the order of arguments and above dependencies

# demo with cylinder
.PHONY: demo
demo: demo/slices.dat
demo/slices.dat: bin/mesh_slicer demo/cylinder.off demo/centreline.dat
	./$^ $@
