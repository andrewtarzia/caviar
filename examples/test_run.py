import stk
import numpy as np
from caviar.cavity_identification import *


class Molecule:

    def __init__(self, stk_molecule):
        self._stk_molecule = stk_molecule

    def select_by_element_string(self, element_string):
        atom_ids = [
            i.get_id() for i in self._stk_molecule.get_atoms()
            if i.__class__.__name__.lower() == element_string.lower()
        ]
        return SelectedMolecule(
            atoms=tuple(self._stk_molecule.get_atoms(atom_ids)),
            position_matrix=np.asarray(tuple(
                self._stk_molecule.get_atomic_positions(atom_ids)
            )),
        )

    def get_num_atoms(self):
        return self._stk_molecule.get_num_atoms()

class SelectedMolecule:

    def __init__(self, atoms, position_matrix):
        self._atoms = atoms
        self._position_matrix = position_matrix

    def get_num_atoms(self):
        return len(self._atoms)

    def get_coords(self):
        return self._position_matrix


class Cavities:

    def __init__(self, cavities, cavities_info, convex_coords):

        self._cavities = cavities
        self._cavities_info = cavities_info
        self._convex_coords = convex_coords

    def _collate_cavities(self):

        coordinates = []
        for i, coords in enumerate(self._cavities):
            c_info = self._cavities_info[i]
            for c, info in zip(coords, c_info):
                coordinates.append((c[0], c[1], c[2], i, info))

        return coordinates

    def write_xyz_file(self, path, cavity='all'):
        """
        Write basic `.xyz` file of Molecule to `path`.

        Connectivity is not maintained in this file type!

        """

        if cavity == 'all':
            coordinates = self._collate_cavities()
        else:
            coordinates = [
                (
                    itm[0], itm[1], itm[2], i,
                    self._cavities_info[cavity][i],
                )
                for i, itm in enumerate(self._cavities[cavity])
            ]

        content = [f'{len(coordinates)}\n\n']
        for coord in coordinates:
            x, y, z, cav_num, decomp = coord
            content.append(
                f'C {x:f} {y:f} {z:f} {cav_num}, {decomp}\n'
            )

        with open(path, 'w') as f:
            f.write(''.join(content))

    def write_convex_xyz_file(self, path):
        """
        Write basic `.xyz` file of Molecule to `path`.

        Connectivity is not maintained in this file type!

        """

        coordinates = []
        for i in self._convex_coords:
            for coord in i:
                coordinates.append(coord)

        content = [f'{len(coordinates)}\n\n']
        for coord in coordinates:
            x, y, z = coord
            content.append(
                f'C {x:f} {y:f} {z:f}\n'
            )

        with open(path, 'w') as f:
            f.write(''.join(content))


def cavity_detection(name, molecule):
    molecule.write(f'{name}_mol.xyz')
    selection_coords = molecule.get_position_matrix()
    selection_molecule = Molecule(molecule)

    # Build a grid around the protein
    grid, grid_shape, grid_min = build_grid(
        selection_coords,
        boxmargin=4.0,
        gridspace=1.0,
    )

    # Main wrapper for cavity identification. A bit of cavity filtering
    # is done here too but it's all geometry based. More will be done
    # afterwards with pharmacophores
    early_cavities, cavities_info, grid_decomposition, array_cavs_coords_hull = (
        wrapper_early_cav_identif(
            name,
            grid,
            grid_min,
            grid_shape,
            selection_molecule,
            selection_coords,
            # size_probe=1.0,
            # maxdistance=20.0,
            radius=3,
            # min_points=2,
            trim_score=0,
            min_burial=7,
            radius_cube = 5,
            # radius_cube_enc = 5,
            min_burial_enc = 2,
        )
    )
    return
    try:
        early_cavities[0]
    except:
        print(
            "CAVIAR does not detect any cavity in with this set of "
            "parameters.\n"
        )

    print(f'there are {len(early_cavities)} cavitites.')

    cavities = Cavities(early_cavities, cavities_info, array_cavs_coords_hull)
    cavities.write_xyz_file('cavity_full.xyz')
    cavities.write_convex_xyz_file('convex.xyz')


def list_of_mols():

    return {
        'cage1': stk.BuildingBlock.init_from_file('cage1.mol'),
        'cage2': stk.BuildingBlock.init_from_file('cage2.mol'),
        'cat': stk.BuildingBlock.init_from_file('cat.mol'),
        'cc3': stk.BuildingBlock.init_from_file('cc3.mol'),
        'complex': stk.BuildingBlock.init_from_file('complex.mol'),
        'moc': stk.BuildingBlock.init_from_file('moc.mol'),
        'big': stk.BuildingBlock.init_from_file('big.mol'),
    }

def main():
    list_of_structures = list_of_mols()

    for s in list_of_structures:
        name = s
        structure = list_of_structures[s]
        cavity_detection(name, structure)


if __name__ == "__main__":
    main()