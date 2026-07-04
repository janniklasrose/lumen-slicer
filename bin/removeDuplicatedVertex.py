#!/usr/bin/env python3
import argparse

import pymeshlab


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a mesh after removing duplicated vertices with PyMeshLab."
    )
    parser.add_argument("input", help="input mesh path")
    parser.add_argument("output", help="output mesh path")
    return parser.parse_args()


def main():
    args = parse_args()

    mesh_set = pymeshlab.MeshSet()
    mesh_set.load_new_mesh(args.input)
    mesh_set.meshing_remove_duplicate_vertices()
    mesh_set.save_current_mesh(args.output, save_vertex_normal=False)


if __name__ == "__main__":
    main()
