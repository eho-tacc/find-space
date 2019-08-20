# !/usr/bin/env python

# Import dependencies
import os
import time
import sys
import numpy as np
import argparse

# ProtPy dependencies
import protpy
from protpy import geom
from protpy.void import VStruct


def set_defaults():
    """Sets default parameters. These values are overwritten by command
    line arguments or JSON file."""

    helpdocs = {}

    helpdocs['work_dir'] = """Path to working directory. FindSpace
    should automatically retrieve this, but I'm keeping it in User
    Defined Parameters in case one wanted more manual control. Ensures
    that all file paths reference from find_space.py, not from the caller's
    working directory."""
    work_dir = "./"

    helpdocs['scan_dir'] = """Path to scan subdirectory
    (need to end in '/'). FindSpace will download all .pdb files it
    screens to scan_dir."""
    scan_dir = "./scan/"

    helpdocs['out_dir'] = """All found void spaces are stored in directories
    named by PDB code in out_dir."""
    out_dir = "./out/"

    helpdocs['pdb_screen_list'] = """File path to pdb screen list, the
    file from which to read PDB codes to be screened."""
    pdb_screen_list = "./pdb_screen_list.txt"

    helpdocs['pdb_list_from'] = """Determines where to get PDB screen
    codes from. Tab-delim will source PDB codes from the first column
    of a tab delimited .txt file. Plain .txt will source
    from a txt file containing a list of PDB codes
    separated by newlines. RCSB will look up the
    current list of all PDB codes in the RCSB and sort them
    alphabetically."""
    # pdb_list_from = "tab_delim"
    pdb_list_from = "plain_txt"
    # pdb_list_from = "rcsb"

    helpdocs['log_fp'] = """File path to log file. This is where FindSpace writes
    PDB codes of structures it has already searched, all
    caps, separated by newlines"""
    log_fp = "./log/fs_log.txt"

    helpdocs['batch_max'] = "Maximum number of structures to scan per run"
    batch_max = 10

    helpdocs['ignore_res'] = """Residues to ignore when computing steric
    hindrances, Voronoi diagram, and convex hull. FindSpace will not
    even bother reading these residues from the PDB file. List of
    residue names separated by commas."""
    ignore_res = 'HOH'

    helpdocs['min_vox_rad'] = """Minimum allowed radius
    (in Angstroms) of pseudoatom voxels."""
    min_vox_rad = 0.5

    helpdocs['max_vox_rad'] = """Maximum allowed radius
    (in Angstroms) of pseudoatom voxels."""
    max_vox_rad = 3.

    helpdocs['bl_mult'] = """Multiplier for minimum allowed voxel size.
    Removes pseudoatom voxels that are within bl_mult * VDW
    radius of the nearest atom in the structure."""
    bl_mult = 1.05

    helpdocs['min_voxel_d'] = """Defines the minimum allowed distance
    between two voxels when FindSpace is making voxels only at Voronoi
    vertices. Defers to a minimum voxel distance of res when voxellating
    along Voronoi ridges. In Angstroms."""
    min_voxel_d = 1.5

    helpdocs['iter_max'] = """FindSpace will cut down on the number of
    pseudo Atom voxels created by iteratively removing redundant voxels
    and removing sterically interfering voxels. This sets a hard
    ceiling for the number of iterations allowed."""
    iter_max = 10

    helpdocs['red_idx'] = """Redundancy index. To remove redundant voxels,
    FindSpace will cycle through all existing voxels and find other
    voxels that are within red_idx * sum of voxel radii in Angstroms.
    It groups these voxels and replaces them with a voxel
    at the arithmetic mean."""
    red_idx = 0.1

    helpdocs['meshvox_res'] = """Meshvoxel resolution. FindSpace renders
    void spaces as meshvoxel points on the surface, and creates a 3D mesh
    as a representation. meshvox_res determines the resolution
    of this mesh representation. Note that higher values,
    especially above 10, will generate more accurate meshes
    but will be more computationally intensive. More accurate
    meshes will give more accurate values for void space volume."""
    meshvox_res = 12.

    helpdocs['env_probe'] = """Environment probe size. FindSpace looks for
    residues in the vicinity of each continuous void space. It will
    probe within env_probe Angstroms of the void space."""
    env_probe = 5.0

    helpdocs['min_void_vol'] = """Sets minimum void volume at which
    FindSpace will report a void space. Set to a value just below the
    estimated VDW volume of your ligand(s) of interest."""
    min_void_vol = 100.

    helpdocs['vdw'] = """Define van der Waals radii for relevant elements,
    in Angstroms. If an atom's VDW is not listed here, will
    default to 2 Angstroms."""
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
    return {
        "work_dir": work_dir,
        "scan_dir": scan_dir,
        "out_dir": out_dir,
        "pdb_screen_list": pdb_screen_list,
        "log_fp": log_fp,
        "pdb_list_from": pdb_list_from,
        "batch_max": batch_max,
        "ignore_res": ignore_res,
        "min_vox_rad": min_vox_rad,
        "max_vox_rad": max_vox_rad,
        "bl_mult": bl_mult,
        "min_voxel_d": min_voxel_d,
        "iter_max": iter_max,
        "red_idx": red_idx,
        "meshvox_res": meshvox_res,
        "env_probe": env_probe,
        "min_void_vol": min_void_vol,
        "vdw": vdw,
        "helpdocs": helpdocs
    }

# end set_defaults


def import_ud_params(defaults):
    """Import command line args to ud_params."""

    p = argparse.ArgumentParser()
    h = defaults['helpdocs']

    p.add_argument("-w", "--work-dir",
                   type=str,
                   default=defaults['work_dir'],
                   help=h['work_dir'])
    p.add_argument("-s", "--scan-dir",
                   help=h['scan_dir'],
                   type=str,
                   default=defaults['scan_dir'])
    p.add_argument("-o", "--out-dir",
                   help=h['out_dir'],
                   type=str,
                   default=defaults['out_dir'])
    p.add_argument("-p", "--pdb-screen-list",
                   help=h['pdb_screen_list'],
                   type=str,
                   default=defaults['pdb_screen_list'])
    p.add_argument("-l", "--log-fp",
                   help=h['log_fp'],
                   type=str,
                   default=defaults['log_fp'])
    p.add_argument("--pdb-list-from",
                   help=h['pdb_list_from'],
                   type=str,
                   choices=['plain_txt', 'tab_delim', 'rcsb'],
                   default='plain_txt')
    p.add_argument("--batch-max",
                   help=h['batch_max'],
                   type=int,
                   default=defaults['batch_max'])
    p.add_argument("--ignore-res",
                   type=str,
                   help=h['ignore_res'],
                   default=defaults['ignore_res'])
    p.add_argument("--min-vox-rad",
                   type=float,
                   help=h['min_vox_rad'],
                   default=defaults['min_vox_rad'])
    p.add_argument("--max-vox-rad",
                   type=float,
                   help=h['max_vox_rad'],
                   default=defaults['max_vox_rad'])
    p.add_argument("--bl-mult",
                   type=float,
                   help=h['bl_mult'],
                   default=defaults['bl_mult'])
    p.add_argument("--min-voxel-d",
                   type=float,
                   help=h['min_voxel_d'],
                   default=defaults['min_voxel_d'])
    p.add_argument("--iter-max",
                   type=int,
                   help=h['iter_max'],
                   default=defaults['iter_max'])
    p.add_argument("--red-idx",
                   type=int,
                   help=h['red_idx'],
                   default=defaults['red_idx'])
    p.add_argument("--meshvox-res",
                   type=int,
                   help=h['meshvox_res'],
                   default=defaults['meshvox_res'])
    p.add_argument("--env-probe",
                   type=int,
                   help=h['env_probe'],
                   default=defaults['env_probe'])
    p.add_argument("--min-void-vol",
                   type=int,
                   help=h['min_void_vol'],
                   default=defaults['min_void_vol'])
    p.add_argument("--vdw",
                   type=str,
                   help="""VDW radii definitions. If you want to change,
                   please define from source (find_space.py line 122)""",
                   default="",
                   dest=None)
    p.add_argument("--print-params",
                   help="Print parameters to console before running",
                   action="store_true")

    args = p.parse_args()
    # TODO: overwrite CLargs with JSON if provided

    # ignore_res should be a list of str
    args.ignore_res = [x.strip().upper() for x in args.ignore_res.split(',')]
    # only use default value for vdw
    args.vdw = defaults['vdw']

    if args.print_params:
        from pprint import PrettyPrinter
        pp = PrettyPrinter()
        print("Runtime parameters:\n")
        pp.pprint(vars(args))
    return vars(args)

# end function import_params


def main_handler(ud_params):
    """
    Main handling function for find_space.py
    """

    batch_max = ud_params['batch_max']
    scan_dir = ud_params['scan_dir']
    out_dir = ud_params['out_dir']
    log_fp = ud_params['log_fp']

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

        """
        ********** LOAD PDB FILE ***********
        """

        # Create Struct object and store in temporary s placeholder
        s = VStruct(pdb, ud_params)

        # DISABLED: Load object into the struct_dict dictionary.
        # Optional step in case comparison between structures will ever
        # be necessary.
        # struct_dict[pdb] = s

        # Downloads PDB file from RCSB to scan_dir. If this fails, delete
        # Struct object and move to next PDB
        if not s.pdbDL(scan_dir):
            os.remove(scan_dir + pdb + ".pdb")  # Remove file from scan_dir
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

            """
            ********** FIND VOID SPACES ***********
            """

            # Find pseudoatom voxels without steric interference and stitch
            # overlapping voxels together
            voids = s.voxelate_voronoi(0,atoms=prot_atoms)

            # Get residues around each void space
            envs = s.void_env(voids,atoms=prot_atoms,byres=True)

            # *******EDITED FOR DEV*******
            """
            ********** SORT VOID SPACES ***********
            """

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

            """
            ********** WRITE VOIDS TO PDB ***********
            """

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
                PDBscribe.write_pdb(void, quiet=True, extra_remark=remark)

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
if __name__ in ("__main__"):
    defaults = set_defaults()
    ud_params = import_ud_params(defaults)
    main_handler(ud_params)
# end if
