import stk

bb1 = stk.BuildingBlock(
    smiles='C1=CC(=CC=C1C2=CC=C(C=C2)N)N',
    functional_groups=[stk.PrimaryAminoFactory()],
)
bb2 = stk.BuildingBlock(
    smiles='C1=C(C=C(C=C1C=O)C=O)C=O',
    functional_groups=[stk.AldehydeFactory()],
)
cage = stk.ConstructedMolecule(
    topology_graph=stk.cage.FourPlusSix(
        building_blocks=(bb1, bb2),
        optimizer=stk.MCHammer(),
    ),
)
cage.write('testmol.mol')