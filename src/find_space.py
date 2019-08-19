# !/usr/bin/env

# Import dependencies
import os, time, sys
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

'''
************** END USER-DEFINED PARAMETERS ***************




**********************************************************
********************** CHANGE LOG ************************
**********************************************************

12/19/2017 - Started development of Find Space, a script designed
    to find void space in protein structures, especially around
    ligands. One of the things I'd like to have while developing
    this is a way to visualize how I'm mapping void space; as
    such, I'm starting by adding a function that can write a
    Struct object to a .pdb file. I can then visualize in
    MacPymol. The next step is to make a set of functions that
    can create pseudoatoms in the structure.
12/21/2017 - Writing a set of functions that can create pseudo-
    atoms in a Struct object. Function voxelate creates pseudoatom
    Atom objects to fill a cube around a given point.
1/10/2018 - I've been flaky about updating the change log. Since
    last time, I started working on the big function of this script,
    Struct.find_space, which is responsible for most of the leg
    work of the script. It analyzes pseudoatom "voxels" created by
    the voxelate function and removes pseudoatoms that are within
    the VDW radius of the nearest atom in the structure (not
    including other pseudoatoms of course). It expands the
    pseudoatoms to the maximum possible size allowed without any
    steric hindrances, and stores this maximum radius arbitrarily
    in the B-factor column of the PDB file. The partner script
    view_space.py translates these B-factors store in the PDB
    and visualizes them as sphere size in the Pymol GUI. Another
    note: Pymol does not recognize atom IDs if they have non-numeric
    characters ("p" prefix will not work). I changed the voxelate
    function to set pseudoatom IDs to the lowest available atom ID.
1/11/2018 - I need a way to exclude pseudoatoms that are in the
    solvated area around the protein structure. The best way to do
    this I think is to construct a convex hull around all the atoms,
    and use that 3D polygon as the border between acceptable and
    excluded pseudoatoms.
1/22/2018 - More or less finished work on voxelate_voronoi method.
    Many updates in the last week:
        -Translates Atom objects into numpy arrays for faster processing.
        This significantly sped up finding each voxel's nearest
        neighbor.
        -Nearest neighbor now uses numpy's KDTree rather than brute force
        -Only Voronoi ridges where both vertices are within the Delaunay
        hull are drawn or considered.
        -User-friendly progress bars!
    Next step is a stitching function that creates sets of pseudoatoms
    that define a void space.
2/27/2018 - Been pretty bad about updating the changelog. I've done
    quite a bit since the last entry.

    1.) voxelate_voronoi often made so many voxel atoms that the .pdb
    file was unopenable. I did a few things to alleviate this:
        -Option in voxelate_line to only voxelate the end points.
        This effectively allows each Voronoi ridge to be represented
        in two voxels or less.
        -Iterate removal of overlapping/redundant voxels:
        Struct.filter_redundant iteratively removes voxels whose centers
        are within red_idx*voxel radius of each other. It replaces these
        two voxels with a new one at the arithmetic mean and recalculates
        its maximum allowed size. This process is repeated a number of
        times as desired by user.
    2.) Added a stitching function. Given a set of voxels, it finds groups
    of them that overlap and assembles them into a list of dicts formatted
    similar to all_atoms [{"id1": obj1,"id2": obj2},{"id3": obj3}].
    3.) At this point, I am trying to find a good way to characterize shape
    and size of a void space (and eventually electrostatics). I started by
    introducing a new type of point, minivoxels, AKA meshvoxels, that are
    spread on the surface of a void space shape. I spent a long time
    figuring out how to calculate a concave hull for the void space,
    only to find out that the concave hull is not at all a good model
    for its volume. I am now working on representing the void space as
    a simpler mesh with triangular facets. This involves several steps:
        -Generate meshvoxels on the outside of the void space
        -Remove meshvoxels (abbr. mvoxs) that lie within, instead of on
        the outside of, voxels representing the void space.
        -Calculate signed volume of the tetrahedron formed by each
        triangular facet and the origin (0,0,0). Apparently, if you
        sum all the signed volumes of these tetrahedrons, you get the
        volume inside the mesh.
3/2/2018 - This past week, I've been working on the Mesh class, which
    is capable of representing irregular spherical shapes as a surface
    mesh. It is capable of importing a dictionary of voxels as a void
    dict, and can also calculate the approximate volume (usually within
    ~5% accuracy) inside the mesh by summing the signed volumes of
    tetrahedrons formed by each simplex and the origin. Today, I will
    try implementing the Mesh class into the find_space script to see
    if it can accurately represent a void space in a PDB file.
3/5/2018 - Successfully implemented Mesh class into FindSpace. I tweaked
    Struct.void_vol so that it will visualize the Mesh as bonded H atoms
    representing meshvoxels. This is the first big milestone for FindSpace;
    it is capable of finding and characterizing void spaces in a molecule
    given a PDB file. Now, it's time to start thinking about future
    directions: what do we do with this information now? So far, we know
    the shape and size of a void space, and can tweak multiple parameters
    in FindSpace to get more accurate, albeit more computationally intensive,
    solutions. These are my thoughts on future directions:
        -It would be nice to have FindSpace export void spaces into more
        permanent storage, so that repeating computationally onorous
        tasks on the whole structure can be avoided. I believe the final
        step of the FindSpace script will be to:
            ~Export the void dict to a separate PDB file entitled
            {PDB}_void. Voxels will be stored as ATOM entries, where
            B-factor will be the radius of the voxel. Each void space
            will have a unique chain designation.
            ~Characterize the environment around the void space. Find
            residues (possibly salts too) in the vicinity of each void
            space, and include this in the _void.pdb file as well.
            What you essentially end up with is a simplified PDB file
            with void spaces included and only relevant, nearby residues
            included.
        -I am able to characterize shape and volume of void spaces. I do not,
        however, plan to use the Mesh class in findspace.py, but rather
        allow the user to perform computationally onorous tasks
        associated with mesh generation in a separate script to be used
        after findspace.py processing. This second script, dedicated to
        characterizing void spaces, will be called ThinkSpace, and will be
        capable of:
            ~Generating Mesh for each void space
            ~Calculating void space volume
            ~Calculate other properties at each meshvoxel, such as
            electrostatic potential
3/7/2018 - Added Struct.void_env function, which identifies residues in
    the vicinity of each void space. It returns an all_atoms like dict
    of Atom objects that are near the void. I implemented the function
    into find_space.py.  Also made a number of small changes throughout
    getfuncy.py and find_space.py:
        1.) Changed all instances of the word "pixel" with "voxel."
        Since what I'm actually making are voxels, not pixels. If you
        see weird variable names like "pnear," it's probably a remnant
        that I didn't catch, and stands for "pixel near."
        2.) Removed reliance on self.all_atoms for all functions. All
        functions that used to use self.all_atoms to reference protein
        atoms (as opposed to voxels) now accept the optional Set param,
        which defaults to self.all_atoms if undefined for the purposes
        of compatibility. Set allows the user to pass in whatever dict
        of Atom objects they like, without having to actually change
        self.all_atoms. The Scribe object has also been edited to write
        a Set dict to PDB file instead of all_atoms from a Structure
        object arg.
        3.) Removed/commented out almost all "DEV" commands, which were
        present to help me debug.
    I tested the script on whole protein PDBs and it does run smoothly,
    albeit slowly. There are still a few problems:
        1.) FindSpace produces more void spaces than are relevant for
        the purposes of most users. I should add a minimum volume of
        void space that must be met before FindSpace reports it. This
        can be done using the Mesh object. This route may be computation-
        -ally heavy, which could possibly be remedied by exiting the mesh
        generation once the void space has reached the minimum volume
        for reporting. Another option is to decrease resolution of the
        Mesh, which will result in an underreported volume, or to estimate
        volume via a sphere overlap algorithm. That is, some areas will
        be subtrated multiple times as areas of overlap, which should
        result in an underestimated void space volume.
        2.) Running the script on large proteins is very slow. The culprits
        are Struct.filter_sterics and Struct.filter_redundant, since they
        are looping over 100,000 times in Python. Consider switching to
        a Numpy-reliant algorithm instead, which should drastically cut down
        on processing time.
3/30/2018 - Main premise of FindSpace is acheived.  Made a new user defined
    parameter that sets a hard minimum on the volume of a void space that
    must be reached before it is reported in the out dir.  I've spent the
    last few days learning about Cython and how to compile it so that I can
    use PyPoint PotentialSolver (https://github.com/SINGROUP/Potential_solver),
    a Cython-based class that solves the Poisson-Boltzmann equation to
    generate an electrostatic potential map.  Before I try to implement this,
    however, I want to get the core functions of FindSpace and getfuncy
    optimized, organized, and well documented, which will take some time.
    Specifically, I am going to restructure the getfuncy module into a
    package hierarchy, with core.py and void.py modules, containing atom-
    related and void space related functions, resp. core.py is imported
    upon import protpy command, and contains essential, universal, non-
    geometry related functions. Submodules such as void and xbond must be
    imported explicitly.
4/3/2018 - Overhauled the file structure of what was formerly known as the
    getfuncy module. I've rearranged this single module into a proto-package
    file hierarchy under the name protpy. It has three submodules:
        -core: Core functions that are dependencies for other submodules.
               These include the Atom, Interaction, and Struct classes,
               housekeeping functions, and basic variable manipulation.
        -geom: Geometric functions. Include Mesh and ConcaveHull classes.
               So far, has no protpy dependencies.
        -void: Space finding functions are now bundled into a Struct
               subclass VStruct, where functions relating to finding and
               manipulating void space are now included as methods.
    I am working on re-writing the VStruct class to optimize performance
    of space finding funcitons. VStruct.filter_sterics has been rewritten
    with heavy Numpy dependence, and is ~5 times faster. I am now working
    on implementing a Cython framework into the protpy module, which will
    (hopefully) become the engine for non-numpy heavy void handling functions.
    I also noticed that relative file paths reference from the caller's
    working directory, not from the location of find_space.py, causing
    several errors when the user runs the script from outside the find_space
    parent dir. I hope to fix this quickly today.
4/6/2018 - Re-wrote filter_steric, filter_red, and stitch_voxels functions in
    Cython. Voxel generation, filtering, and stitching is now literally
    thousands of times faster.  The geom.Mesh class is still very slow, and is,
    I believe, the last remaining step in the core FindSpace process that is
    annoyingly slow. One confusing step is how to determine handedness of a
    triangular simplex in Cython. Before, I calculated the Delaunay hull and
    used in_hull function to determine if the normal was inside the mesh. This
    would not work in Cython. Instead, consider determining if normal vector is
    pointing right way by calculating two possible normals and asking which one
    is closer to the center of the voxel.
4/9/2018 - Re-wrote a portion of Mesh.import_voxels in Cython. The only really
    slow part of this function was finding simplices of the mesh where one of
    the vertices laid inside the mesh. Not really sure why this process was
    so slow in Python, but it is now hundreds of times faster with the Cython
    implementation, remove_inner_simps. I changed the red_idx UD param to
    mark voxels as redundant if they are within red_idx times the sum of the
    voxel radii.  FindSpace now runs pretty fast, processing the relatively
    large structure 1GSK in about 2 minutes, most of which was used voxellating
    along Voronoi ridges (still a Python implementation). I am satisfied with
    how fast the script is currently.  I will update docs on the ProtPy library
    and maybe add some other user friendly functionality, but otherwise let
    this script be for a while and start working on mapping electrostatics in
    the output PDB files that FindSpace produces.
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

    # Just runs the script once
    main_handler()

# end if
