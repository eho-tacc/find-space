# !/usr/bin/env python

# Import dependencies
import os
import time
import sys
import numpy as np

# ProtPy dependencies
import protpy
from protpy import geom
from protpy.void import VStruct
#from electrostatic_potential_solver import ElectrostaticPotentialSolver

'''
*********************************************************
**************** USER-DEFINED PARAMETERS ****************
*********************************************************
'''

# Defines the base directory. FindSpace should automatically
# retrieve this, but I'm keeping it in User Defined Parameters
# in case one wanted more manual control. Ensures that all
# file paths reference from find_space.py, not from the caller's
# working directory.
base_dir = "./" + os.path.join(os.path.dirname(__file__), '')

# Location of out and scan subdirectories (need to end in '/').
# FindSpace will download all .pdb files it screens to scan_dir.
# All found void spaces are stored in directories named by PDB code
# in out_dir.
scan_dir = base_dir + "scan/"
out_dir = base_dir + "out/"

# File path of pdb screen list, the file from which to
# read PDB codes to be screened.
pdb_screen_list = base_dir + "pdb_screen_list.txt"

# Determines where to get PDB screen codes from.
# Uncomment only one of the following lines. Tab-
# delim will source PDB codes from the first column
# of a tab delimited .txt file. Plain .txt will source
# from a txt file containing a list of PDB codes
# separated by newlines (\n). RCSB will look up the
# current list of all PDB codes and sort them
# alphabetically.
# pdb_list_from = "tab_delim"
pdb_list_from = "plain_txt"
# pdb_list_from = "rcsb"

# File path to log file. This is where FindSpace writes
# PDB codes of structures it has already searched, all
# caps, separated by newlines (\n)
log_fp = base_dir + "log/fs_log.txt"

# Define the maximum number of structures to scan per run
batch_max = 10

# Residues to ignore when computing steric hindrances,
# Voronoi diagram, and convex hull. FindSpace will not
# even bother reading these residues from the PDB file.
ignore_res = ['HOH']

# Minimum and maximum allowed radius (in Angstroms) of pseudoatom voxels.
vox_rad_lims = (0.5,3.)

# Multiplier for minimum allowed voxel size. Excludes
# pseudoatom voxels that are within bl_mult * VDW
# radius of the nearest atom in the structure.
bl_mult = 1.05

# Defines the minimum allowed distance between two voxels
# when FindSpace is making voxels only at Voronoi vertices.
# Defers to a minimum voxel distance of res when voxellating
# along Voronoi ridges. In Angstroms.
min_voxel_d = 1.5

# FindSpace will cut down on the number of pseudo Atom voxels
# created by iteratively removing redundant voxels and
# removing sterically interfering voxels. This sets a hard
# ceiling for the number of iterations allowed.
iter_max = 10

# Redundancy index. To remove redundant voxels, FindSpace
# will cycle through all existing voxels and find other
# voxels that are within red_idx * sum of voxel radii in Angstroms.
# It groups these voxels and replaces them with a voxel
# at the arithmetic mean.
red_idx = 0.1

# Meshvoxel resolution. FindSpace renders void spaces as
# meshvoxel points on the surface, and creates a 3D mesh
# as a representation. meshvox_res determines the resolution
# of this mesh representation. Note that higher values,
# especially above 10, will generate more accurate meshes
# but will be more computationally intensive. More accurate
# meshes will give more accurate values for void space volume.
meshvox_res = 12.

# Environment probe. FindSpace looks for residues in the vicinity
# of each continuous void space. It will probe within env_probe
# Angstroms of the void space.
env_probe = 5.0

# Sets minimum void volume at which FindSpace will report
# a void space. Set to a value just below the estimated
# VDW volume of your ligand(s) of interest.
min_void_vol = 100.

# Define van der Waals radii for relevant elements, in Angstroms.
# If an atom's VDW is not listed here, X-finder will
# default to 2 Angstroms.
vdw = {
    'H': 1.2,
    'F': 1.5,
    'CL': 1.9,
    'Br': 1.97,
    'I': 2.15,
    'O': 1.58,
    'S': 1.85,
    'N': 1.64,
    'C': 1.77,
    'default': 2
}

# Package all user-defined parameters into a dictionary
ud_params = {
    "vdw": vdw,
    "pdb_screen_list": pdb_screen_list,
    "log_fp": log_fp,
    "pdb_list_from": pdb_list_from,
    "batch_max": batch_max,
    "ignore_res": ignore_res,
    "vox_rad_lims": vox_rad_lims,
    "bl_mult": bl_mult,
    "min_voxel_d": min_voxel_d,
    "iter_max": iter_max,
    "red_idx": red_idx,
    "meshvox_res": meshvox_res,
    "env_probe": env_probe,
    "min_void_vol": min_void_vol
}


def import_params():
    '''
    Import command line args to ud_params.
    '''
    pass
# end function import_params


'''
************** END USER-DEFINED PARAMETERS ***************
'''


def main_handler():
    '''
    Main handling function for find_space.py
    '''

    # Define dictionary of PDB Struct objects
    struct_dict = {}

    # Define list of PDB IDs that have been searched this round
    searched_this_round = []

    # Defines list of PDB codes to screen
    pdbs_to_screen = []

    # Retrieves the list of PDB codes to screen, either from the
    # pdb_screen_list file or from the RCSB PDB datbase.
    pdbs_to_screen = protpy.get_screen_codes(ud_params)

    # Exit script if all PDBs have already scanned
    if not pdbs_to_screen:
        print("*  main_handler - All PDBs listed in ./log/search_log.txt" +
              " have been screened.  If you wish to screen these again," +
              " please remove the appropriate PDB codes from the" +
              " ./log/search_log.txt file.")
        return
    # end if

    # Loop through download list and look for halogens in each PDB file
    # one at a time. Limit to batch_max structures per cycle
    for pdb in pdbs_to_screen[0:batch_max]:

        # Instantiate timer
        pdb_timer = protpy.Timer()

        # Create Scribe object in charge of logging searched PDBs
        # in the search_log.txt file.
        Log = protpy.Scribe("./log/search_log.txt",pdb + "\n")

        '''
        ********** LOAD PDB FILE ***********
        '''

        # Create Struct object and store in temporary s placeholder
        s = VStruct(pdb,ud_params)

        # DISABLED: Load object into the struct_dict dictionary.
        # Optional step in case comparison between structures will ever
        # be necessary.
        # struct_dict[pdb] = s

        # Downloads PDB file from RCSB to scan_dir. If this fails, delete
        # Struct object and move to next PDB
        if not s.pdbDL(scan_dir):
            os.remove(scan_dir + pdb + ".pdb") # Remove file from scan_dir
            Log.write(append=True)
            del s
            continue
        # end if

        # Reads PDB file, instantiates Atom objects,
        # and logs bonds between Atoms
        prot_atoms = s.read_pdb(scan_dir + pdb + ".pdb",ud_params)

        # Wraps the computationally intensive stuff in a try statement
        # with a KeyboardInterrupt exception in case the
        # computation gets really long and you want to do
        # something else with your life besides watching
        # Terminal.
        try:

            '''
            ********** FIND VOID SPACES ***********
            '''

            # Find pseudoatom voxels without steric interference and stitch
            # overlapping voxels together
            voids = s.voxelate_voronoi(0,atoms=prot_atoms)

            # Get residues around each void space
            envs = s.void_env(voids,atoms=prot_atoms,byres=True)

            # *******EDITED FOR DEV*******
            '''
            ********** SORT VOID SPACES ***********
            '''

            # Calculate the volumes of each void space
            vols = s.void_vol(voids)

            # Convert voids, envs, and vols to ndarray
            voids = np.array(voids)
            vols = np.array(vols)
            envs = np.array(envs)

            # Indices of void spaces with vols > min_void_vol
            big_enough = np.where(vols>ud_params['min_void_vol'])[0]

            # Filters out void spaces in vols, envs, and voids
            # that are not big enough
            voids = voids[big_enough]
            vols = vols[big_enough]
            envs = envs[big_enough]

            '''
            ********** WRITE VOIDS TO PDB ***********
            '''

            # Double check to make sure there are the same
            # number of elements in voids and envs
            if len(voids) != len(envs):
                print("   main_handler - Error finding " + \
                      "protein atoms around void spaces.")
                continue
            # end if

            # User friendly
            print("   main_handler - Writing void spaces " + \
                  "to PDB files...")
            write_bar = protpy.ProgressBar(len(voids))

            # Loop through all voids and associated envs
            for void,env,ct in \
                zip(voids,envs,range(len(voids)+1)):

                # Instantiate Scribe object
                # responsible for writing .pdb file
                PDBscribe = protpy.Scribe(
                    out_dir + pdb + "/" + pdb + \
                    "_void" + str(ct) + ".pdb"
                    )

                # Add env Atom dict to void dict
                void.update(env)

                # Extra remark in PDB file
                remark = pdb + "_void" + str(ct) + " volume=" + \
                str('{0:.1f}'.format(vols[ct])) + " Ang^3, " + \
                str(len(void)) + " voxels, " + str(len(env)) + \
                " environment atoms"

                # Import void Atom dict into
                # PDBscribe and format in .pdb file format
                PDBscribe.write_pdb(void,quiet=True,
                                    extra_remark=remark)

                # Write PDBscribe text to PDB file
                PDBscribe.write()

                # Update prog bar
                write_bar.update()
            # end for

            # Wrap it up: write pdb to log file,
            # delete PDB file, and delete Struct obj
            Log.write(append=True)
            del s

        except KeyboardInterrupt:

            # user got bored of waiting and aborted
            # with ctrl-C
            print(
                "\n*  main_handler - Script aborted by user. " + \
                pdb + " will not be added to search log."
            )
            break

        finally:

            # Print timer info
            print(
                "   main_handler - FindSpace spent " + \
                str('{0:.0f}'.format(pdb_timer.stop())) + \
                " seconds on finding voids in " + pdb + ". "
            )

            # Clean search log
            protpy.clean_search_log(log_fp)

        # end try

    # end for

# end function main_handler

# __name__ == "pymol" if running from PyMol command GUI
if __name__ in ("__main__","pymol"):

    parser = argparse.ArgumentParser()
    import_params(parser)
    main_handler()

# end if
