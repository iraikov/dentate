import itertools

import h5py

h5types_file = 'dentate_h5types.h5'

DG_populations = ["AAC", "BC", "GC", "HC", "HCC", "IS", "MC", "MOPP", "NGFC", "MPP", "LPP", "CA3c", "ConMC"]
DG_IN_populations = ["AAC", "BC", "HC", "HCC", "IS", "MC", "MOPP", "NGFC"]
DG_EXT_populations = ["MPP", "LPP", "CA3c", "ConMC"]

DG_coords_file = "DG_coords_20190717.h5"

DG_GC_coordinate_file  = "DGC_coords_reindex_20190717.h5"
DG_IN_coordinate_file  = "dentate_IN_coords_dist_20190622.h5"
DG_EXT_coordinate_file = "dentate_IN_coords_dist_20190622.h5"


coordinate_files = {
     'AAC':  DG_IN_coordinate_file,
     'BC':   DG_IN_coordinate_file,
     'GC':   DG_GC_coordinate_file,
     'HC':   DG_IN_coordinate_file,
     'HCC':  DG_IN_coordinate_file,
     'IS':   DG_IN_coordinate_file,
     'MC':   DG_IN_coordinate_file,
     'MOPP': DG_IN_coordinate_file,
     'NGFC': DG_IN_coordinate_file,
     'MPP':  DG_EXT_coordinate_file,
     'LPP':  DG_EXT_coordinate_file,
     'CA3c':  DG_EXT_coordinate_file,
     'ConMC':  DG_EXT_coordinate_file
}

coordinate_ns_generated = 'Generated Coordinates'
coordinate_ns_interpolated = 'Interpolated Coordinates'
coordinate_namespaces = {
     'AAC':  coordinate_ns_generated,
     'BC':   coordinate_ns_generated,
     'GC':   coordinate_ns_interpolated,
     'HC':   coordinate_ns_generated,
     'HCC':  coordinate_ns_generated,
     'IS':   coordinate_ns_generated,
     'MC':   coordinate_ns_generated,
     'MOPP': coordinate_ns_generated,
     'NGFC': coordinate_ns_generated,
     'MPP':  coordinate_ns_generated,
     'LPP':  coordinate_ns_generated,
     'CA3c':  coordinate_ns_generated,
     'ConMC':  coordinate_ns_generated
}

## Creates H5Types entries
with h5py.File(DG_coords_file) as f:
    input_file  = h5py.File(h5types_file,'r')
    input_file.copy('/H5Types',f)
    input_file.close()

## Creates coordinates entries
with h5py.File(DG_coords_file) as f:

    grp = f.create_group("Populations")
                
    for p in DG_populations:
        grp.create_group(p)

    for p in DG_populations:
        coords_file = coordinate_files[p]
        coords_ns   = coordinate_namespaces[p]
        group_id    = f.require_group('/Populations/%s' % p)

        with h5py.File(coords_file) as fs:
            fs.copy('/Populations/%s/%s' % (p, coords_ns), group_id, name='Coordinates')
            if ('/Populations/%s/Arc Distances' % p) in fs:
                fs.copy('/Populations/%s/Arc Distances' % p, group_id, name='Arc Distances')
