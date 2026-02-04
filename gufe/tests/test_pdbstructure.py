import pytest
from io import StringIO
from gufe.vendor.pdb_file.pdbstructure import PdbStructure, _parse_atom_index

def test_hex_conect_parsing():
    pdb_snippet = """ATOM  99999  C   LIG A   1       0.000  0.000  0.000  1.00 0.00           C
ATOM  A000F  N   LIG A   2       1.000  0.000  0.000  1.00 0.00           N
ATOM  A000G  O   LIG A   3       0.000  1.000  0.000  1.00 0.00           O
CONECT99999A000FA000G"""
    print(pdb_snippet)

    f = StringIO(pdb_snippet)
    pdb = PdbStructure(f, load_all_models=True)

    # Collect all atom serial numbers, including Maestro-style
    atom_serials = [atom.serial_number for atom in
                    pdb.iter_atoms(use_all_models=True)]

    # There should be 3 atoms
    assert len(atom_serials) == 3

    # All serial numbers should be integers and unique
    for serial in atom_serials:
        assert isinstance(serial, int)
    assert len(set(atom_serials)) == 3

    # Convert the known Maestro-style indices to integers
    a000f_serial = _parse_atom_index("A000F")
    a000g_serial = _parse_atom_index("A000G")
    assert a000f_serial in atom_serials
    assert a000g_serial in atom_serials

    # Check that CONECT records refer to the correct integers
    conects = pdb._current_model.connects
    print('CONECTS', conects)
    assert len(conects) == 1
    central, bonded1, bonded2 = conects[0]
    assert central == 99999
    # The bonded atoms match the converted Maestro-style serials
    assert bonded1 in [a000f_serial, a000g_serial]
    assert bonded2 in [a000f_serial, a000g_serial]
    assert bonded1 != bonded2

    # Optional: test that specific known hex serials were converted correctly
    known_hex_serials = ["A000G", "A000F"]  # replace with actual serials in your file
    for hex_serial in known_hex_serials:
        idx = _parse_atom_index(hex_serial)
        assert isinstance(idx, int)
        assert idx >= 100000  # confirms it went through hex -> int logic
