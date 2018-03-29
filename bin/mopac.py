

def make_input_file(atoms, coordinates, header, charge, parameters):
    """
CHARGE=0 PM6-D3 SINGLET PRECISE NOMM 
g0-0.out

C  -0.00000 1  0.00001 1  0.00000 1
H   0.60979 1  0.30198 1 -0.85268 1
H  -0.29703 1 -1.04321 1 -0.11664 1
H  -0.89050 1  0.62812 1  0.05088 1
H   0.57775 1  0.11307 1  0.91843 1

    """


    out = ""
    coordinate_line = "{:2s}  {:13.10f} 1  {:13.10f} 1  {:13.10f} 1\n"

    # fix charge
    if not charge == 0:
        header = header.replace("charge=0", "charge="+str(charge))

    # energy model
    out += header

    out += "title" + "\n"
    out += "\n"

    # Coordinates
    for atom, coord in zip(atoms, coordinates):
        out += coordinate_line.format(atom, *coord)

    return out

