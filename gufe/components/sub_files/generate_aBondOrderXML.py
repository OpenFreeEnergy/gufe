import xml.etree.ElementTree as etree

out_path = "data/residues.xml"
in_path = "data/residues_orig.xml"

import xml.etree.ElementTree as etree

out_path = "../gufe/components/sub_files/data/residues.xml"
in_path = "../gufe/components/sub_files/data/residues_orig.xml"


exception_bond_keys = {
            # AminoAcids
            ##Backbone
            ('C', 'O'): { "order": 2, "resns": "all"},
            ## Carbonyls in R
            ("CZ", "NH2"): { "order": 2, "resns": ("ARG")},
            ("CG", "OD1"):{ "order": 2, "resns":  ("ASP", "ASN")},
            ("CD", "OE1"):{ "order": 2, "resns": ("GLN", "GlU")},
            ("CD", "OE1"):{ "order": 2, "resns": ("GlU")},
            ("CD", "OE"):{ "order": 2, "resns": ("PCA")},
            
            ## Aromatics:
            ("CD", "CG"):{ "order": 2, "resns": ("HIS")},
            ("CE1", "NE2"):{ "order": 2, "resns": ("HIS")},
            
            ("CG", "CD1"):{ "order": 2, "resns": ("PHE")},
            ("CE1", "CZ"):{ "order": 2, "resns": ("PHE")},
            ("CE2", "CD2"):{ "order": 2, "resns": ("PHE")},
            
            ("CD1", "NE1"):{ "order": 2, "resns": ("TRP")},
            ("CD2", "CE3"):{ "order": 2, "resns": ("TRP")},
            ("CE2", "CZ2"):{ "order": 2, "resns": ("TRP")},
            ("CZ3", "CH2"):{ "order": 2, "resns": ("TRP")},

            ("CG", "CD1"):{ "order": 2, "resns": ("TYR")},
            ("CE1", "CZ"):{ "order": 2, "resns": ("TYR")},
            ("CD2", "CE2"):{ "order": 2, "resns": ("TYR")},

            # NucleicAcids
            ## Phosphates
            ("OP1", "P"):{ "order": 2, "resns": ("U", "G", "A", "C", "DT", "DG", "DC", "DA")},
            
            ## Pyrimidines: Uracil, Thymin and Cytosin
            ("C2", "O2"):{ "order": 2, "resns": ("U", "DT", "C", "DC")},
            ("C5", "C6"):{ "order": 2, "resns": ("U", "DT", "C", "DC")},
            ("C4", "O4"):{ "order": 2, "resns": ("U", "DT")},
            ("C4", "N3"):{ "order": 2, "resns": ("C", "DC")},

            ## Purines: Guanine, Adenine
            ("C2", "N3"):{ "order": 2, "resns": ("G", "DG", "A", "DA")},
            ("C4", "C5"):{ "order": 2, "resns": ("G", "DG", "A", "DA")},
            ("N7", "C8"):{ "order": 2, "resns": ("G", "DG", "A", "DA")},
            ("C6", "O6"):{ "order": 2, "resns": ("G", "DG")},
            ("C6", "N1"):{ "order": 2, "resns": ("A", "DA")},            
             }


tree = etree.parse(in_path)

for residue in tree.getroot().findall('Residue'):
    resn = residue.get("name")
    #print(residue.get("name"))
    for bond in residue.findall("Bond"):
        #print(bond)
        c1 = bond.get("from")
        c2 = bond.get("to")
        bond_atoms=tuple(sorted([c1, c2]))
        
        if(bond_atoms in exception_bond_keys and (not "resn" in exception_bond_keys[bond_atoms] or exception_bond_keys[bond_atoms]["resn"] == resn)):
            bond.set("order", str(exception_bond_keys[bond_atoms]["order"]))
        else:
            bond.set("order", str(1))
tree.write(out_path)