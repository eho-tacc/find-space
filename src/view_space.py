# !/usr/bin/env

'''
**********************************************************
********************** CHANGE LOG ************************
**********************************************************

1/9/2018 - Started working on view_space.py, a Pymol
    exclusive script that serves as a viewer for the results
    of find_space.py. Is able to read PDB files written by
    find_space.py and converts B-factors of pseudoatoms to
    the sphere_scale attribute in Pymol. This allows for
    scaled visualization of empty space in the structure.
'''

def vdwToB(ID,b,vdw_max=3):
    '''
    Sets VDW radius to the B-factor for atom with ID.

    Arguments:
        int ID - ID number of atom to alter
        float b - B factor of the atom
        float vdw_max - maximum shown size of VDW radius
    '''

    # Limit VDW radius so spheres don't explode
    if b > vdw_max:
        b = vdw_max
    # end if

    cmd.alter("id " + str(ID), "vdw=" + str(b))

# end function test_func

def main_handler():
    '''
    Main handling function for view_space.py.
    '''

    cmd.select("pixels","chain X")
    cmd.select("prot","not pixels")
    cmd.iterate("pixels","vdwToB(ID,b)")
    cmd.hide("everything")
    cmd.show("lines","chain X and elem H")
    cmd.show("spheres","chain X")
    cmd.show("lines","prot")
    cmd.show("sticks","chain J")
    cmd.set("sphere_transparency","0.3")
    cmd.rebuild()

# end function main_handler

# __name__ == "pymol" if running from PyMol command GUI
if __name__ in ("pymol"):
    main_handler()
# script was called from an environment other than PyMol
else:
    print("   x_review.py - Please run X-review from " +
          "the PyMol command line.")
# end if
