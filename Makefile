# by default, compile only
default: bin/mesh_slicer

clean:
	$(RM) -f bin/mesh_slicer

# compile the slicing executable
CXXFLAGS = -O2 -std=c++11
bin/mesh_slicer: src/mesh_slicer.cpp
	$(CXX) $(CPPFLAGS) $< -o $@ -lCGAL -lgmp
