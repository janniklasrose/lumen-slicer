# Lumen Slicer

Methods to analyse lumina (tubular geometries) by slicing along their centreline

## Getting started

Create and activate a Python virtual environment:

```sh
python3 -m venv .venv
source .venv/bin/activate
```

Install PyMeshLab:

```sh
python -m pip install pymeshlab
```

Install CGAL:

```sh
brew install cgal
```

Compile the tools:

```sh
make
```

## Tools

- `mesh_slicer mesh.off centreline.dat slices.dat` slices a triangle mesh with the planes in `centreline.dat`.
- `tessellate_slices slices.dat slices.off` converts slice contours into a tessellated OFF mesh.
  - Pass `--refine` to `tessellate_slices` to enable CGAL mesh refinement with default criteria.
- `surface_reconstructor slices.dat mesh.off` reconstructs an OFF triangle mesh from the slice contours written by `mesh_slicer`.

## Demo

Run the cylinder demo to generate `demo/slices.dat`, `demo/slices.off`, and `demo/reconstructed.off`:

```sh
make demo
```
