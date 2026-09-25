"""Measured PDB features for single-chain, standard-residue peptide candidates.

This module does not predict or repair structures. It rejects mismatched or
incomplete inputs instead of manufacturing coordinates or missing features.
"""
import copy

import numpy as np
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.Data.PDBData import protein_letters_3to1

SIDECHAIN = {
    "A": "CB", "C": "CB SG", "D": "CB CG OD1 OD2", "E": "CB CG CD OE1 OE2",
    "F": "CB CG CD1 CD2 CE1 CE2 CZ", "G": "", "H": "CB CG ND1 CD2 CE1 NE2",
    "I": "CB CG1 CG2 CD1", "K": "CB CG CD CE NZ", "L": "CB CG CD1 CD2",
    "M": "CB CG SD CE", "N": "CB CG OD1 ND2", "P": "CB CG CD",
    "Q": "CB CG CD OE1 NE2", "R": "CB CG CD NE CZ NH1 NH2", "S": "CB OG",
    "T": "CB OG1 CG2", "V": "CB CG1 CG2", "W": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "Y": "CB CG CD1 CD2 CE1 CE2 CZ OH",
}


def read_peptide(path, sequence, chain_id, model_index=1):
    structure = PDBParser(PERMISSIVE=False, QUIET=True).get_structure("candidate", str(path))
    models = list(structure)
    if type(model_index) is not int or not 1 <= model_index <= len(models):
        raise ValueError("model_index is 1-based and must select an existing PDB model")
    model = models[model_index - 1]
    if not isinstance(chain_id, str) or chain_id not in model:
        raise ValueError("Explicit existing PDB chain_id required")
    # Context-dependent features must not silently drop other chains or ligands.
    if len(list(model)) != 1:
        raise ValueError("Supply an isolated single-chain peptide PDB; multichain context is unsupported")
    chain = copy.deepcopy(model[chain_id])
    residues = list(chain)
    if any(r.id[0] != " " or r.resname not in protein_letters_3to1 for r in residues):
        raise ValueError("Only standard unmodified residues supported; prepare water/ligand-free peptide PDB")
    actual = "".join(protein_letters_3to1[r.resname] for r in residues)
    if actual != sequence:
        raise ValueError(f"Candidate PDB sequence mismatch: expected {sequence}, got {actual}")
    for residue, aa in zip(residues, sequence):
        if residue.is_disordered():
            raise ValueError("Resolve alternate conformations before importing the structure")
        required = {"N", "CA", "C", "O", *SIDECHAIN[aa].split()}
        missing = required - {a.name for a in residue}
        if missing:
            raise ValueError(f"Missing heavy atoms at {residue.id}: {sorted(missing)}")
        for atom in list(residue):
            if not np.isfinite(atom.coord).all():
                raise ValueError("Nonfinite PDB coordinates")
            if atom.element in ("H", "D"):
                residue.detach_child(atom.id)
            elif atom.element not in ("C", "N", "O", "S"):
                raise ValueError(f"Unexpected atom element {atom.element}")
            elif atom.occupancy is None or atom.occupancy <= 0:
                raise ValueError("Missing or zero heavy-atom occupancy")
    return chain, residues


def features(path, sequence, position, chain_id, model_index=1):
    if type(position) is not int or not 1 <= position <= len(sequence) or sequence[position - 1] != "C":
        raise ValueError("Target must be a valid 1-based Cys")
    chain, residues = read_peptide(path, sequence, chain_id, model_index)
    atoms = list(chain.get_atoms())
    coords = np.array([a.coord for a in atoms], dtype=float)
    target = residues[position - 1]
    sg = np.array(target["SG"].coord, dtype=float)
    sr = ShrakeRupley(probe_radius=1.4, n_points=960)
    sr.compute(chain, level="A")
    sasa = float(target["CB"].sasa + target["SG"].sasa)
    distances = np.linalg.norm(coords - sg, axis=1)
    nearby = int(sum(d <= 5.0 and a is not target["SG"] for a, d in zip(atoms, distances)))
    sulfur = [{"position_1based": j + 1, "distance_A": float(np.linalg.norm(r["SG"].coord - sg))}
              for j, r in enumerate(residues) if r is not target and r.resname == "CYS"]
    sulfur.sort(key=lambda r: r["distance_A"])
    # Geometry diagnostics are broad validity checks, not a structure-energy model.
    backbone = []
    for j, r in enumerate(residues):
        for a, b, low, high in [("N", "CA", 1.2, 1.7), ("CA", "C", 1.2, 1.7), ("C", "O", 1.0, 1.5)]:
            distance = float(np.linalg.norm(r[a].coord - r[b].coord))
            if not low <= distance <= high:
                backbone.append({"position_1based": j + 1, "bond": f"{a}-{b}", "distance_A": distance})
        if j:
            distance = float(np.linalg.norm(residues[j-1]["C"].coord - r["N"].coord))
            if not 1.1 <= distance <= 1.6:
                backbone.append({"position_1based": j + 1, "bond": "previous_C-N", "distance_A": distance})
    for a, b, low, high in [("CA", "CB", 1.2, 1.8), ("CB", "SG", 1.5, 2.1)]:
        distance = float(np.linalg.norm(target[a].coord - target[b].coord))
        if not low <= distance <= high:
            backbone.append({"position_1based": position, "bond": f"{a}-{b}", "distance_A": distance})
    atom_residue = [residues.index(a.parent) for a in atoms]
    clashes = []
    for i, a in enumerate(atoms):
        # Severe overlaps, excluding same-residue bonded atoms and sequential residues.
        ri = atom_residue[i]
        for j in range(i + 1, len(atoms)):
            rj = atom_residue[j]
            if abs(ri - rj) <= 1:
                continue
            d = float(np.linalg.norm(coords[i] - coords[j]))
            if d < 1.2:
                clashes.append({"residue_a": ri + 1, "atom_a": a.name,
                                "residue_b": rj + 1, "atom_b": atoms[j].name, "distance_A": d})
    length = len(sequence)
    result = {"seq": sequence, "xyz": [r["CA"].coord.astype(float).tolist() for r in residues],
              **{k: [0] * length for k in ("is_disulfide", "sasa", "prev_is_pro", "nearby_atoms", "feature_mask")}}
    idx = position - 1
    result["is_disulfide"][idx] = int(any(r["distance_A"] < 2.5 for r in sulfur))
    result["sasa"][idx] = sasa
    result["prev_is_pro"][idx] = int(idx > 0 and sequence[idx - 1] == "P")
    result["nearby_atoms"][idx] = nearby
    result["feature_mask"][idx] = 1
    evidence = {"target_pdb_residue": {"chain": chain_id, "number": target.id[1], "insertion_code": target.id[2]},
                "model_index": model_index, "sasa_sidechain_A2": sasa,
                "sasa_method": "Biopython ShrakeRupley; heavy atoms; CB+SG; probe 1.4 A; 960 points",
                "nearby_atoms_5A": nearby, "nearby_definition": "heavy atoms within 5 A of target SG, excluding target SG only",
                "other_cys_SG_distances": sulfur, "backbone_geometry_violations": backbone,
                "severe_nonlocal_overlaps": clashes, "geometry_passed": not backbone and not clashes,
                "scope": "geometric checks only, not folding accuracy or reaction feasibility",
                "feature_compatibility": "historical training-feature computation not fully documented; distribution shift possible"}
    return result, evidence, chain


def write_chain(chain, destination):
    writer = PDBIO()
    writer.set_structure(chain)
    writer.save(str(destination))
