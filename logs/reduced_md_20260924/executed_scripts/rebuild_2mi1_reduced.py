"""Build reduced 2MI1 and run a fixed short, multi-start explicit-water MD protocol.

Run with OpenMM 8.6.1 and PDBFixer 1.12.0 installed, or use the isolated
tmp/agents/md_runtime directory. No scoring result controls sampling.
"""
import argparse
import io
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
if (ROOT/'tmp/agents/md_runtime').is_dir():
    sys.path.insert(0, str(ROOT/'tmp/agents/md_runtime'))
import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer
from src import core


def no_disulfide(topology, positions):
    """Rebuild topology with no hydrogens or SG-SG bonds; preserve all other bonds."""
    new = app.Topology()
    mapping, coordinates, removed = {}, [], []
    for chain in topology.chains():
        c = new.addChain(chain.id)
        for residue in chain.residues():
            r = new.addResidue(residue.name, c, residue.id, residue.insertionCode)
            for a in residue.atoms():
                if a.element == app.element.hydrogen:
                    continue
                mapping[a] = new.addAtom(a.name, a.element, r, a.id)
                coordinates.append(positions[a.index].value_in_unit(unit.nanometer))
    for a, b in topology.bonds():
        if a not in mapping or b not in mapping:
            continue
        if a.name == b.name == 'SG':
            removed.append([a.residue.id,b.residue.id])
        else:
            new.addBond(mapping[a], mapping[b])
    return app.Modeller(new, unit.Quantity(coordinates, unit.nanometer)), removed


def write_pdb(path, topology, positions):
    with path.open('w') as stream:
        app.PDBFile.writeFile(topology, positions, stream, keepIds=True)


def cap_termini(modeller):
    """Add ACE carbonyl/methyl and NHE primary amide heavy atoms; relax afterwards."""
    old=modeller.topology
    xyz=np.asarray(modeller.positions.value_in_unit(unit.angstrom))
    residues=list(old.residues())
    atom=lambda r,name: next(a for a in r.atoms() if a.name==name)
    normalize=lambda v:v/np.linalg.norm(v)
    n=xyz[atom(residues[0],'N').index];ca=xyz[atom(residues[0],'CA').index]
    c=xyz[atom(residues[0],'C').index]
    v=normalize(ca-n);p=normalize((c-ca)-np.dot(c-ca,v)*v)
    candidates=[n+1.33*(-0.5*v+s*np.sqrt(3)/2*p) for s in (-1,1)]
    capc=max(candidates,key=lambda q:np.min(np.linalg.norm(xyz[2:]-q,axis=1)))
    w=normalize(n-capc);p=normalize((ca-n)-np.dot(ca-n,w)*w)
    capo=capc+1.23*(-0.5*w+np.sqrt(3)/2*p)
    methyl=capc+1.50*(-0.5*w-np.sqrt(3)/2*p)
    last=residues[-1];c=xyz[atom(last,'C').index];ca=xyz[atom(last,'CA').index];o=xyz[atom(last,'O').index]
    amiden=c+1.335*normalize(-normalize(ca-c)-normalize(o-c))
    top=app.Topology();chain=top.addChain('A');coords=[];mapping={}
    ace=top.addResidue('ACE',chain,'0')
    ace_atoms={}
    for name,element,position in [('CH3',app.element.carbon,methyl),('C',app.element.carbon,capc),('O',app.element.oxygen,capo)]:
        ace_atoms[name]=top.addAtom(name,element,ace);coords.append(position)
    top.addBond(ace_atoms['CH3'],ace_atoms['C']);top.addBond(ace_atoms['C'],ace_atoms['O'])
    for r in residues:
        nr=top.addResidue(r.name,chain,r.id)
        for a in r.atoms():
            if a.name=='OXT':continue
            mapping[a]=top.addAtom(a.name,a.element,nr);coords.append(xyz[a.index])
    for a,b in old.bonds():
        if a in mapping and b in mapping:top.addBond(mapping[a],mapping[b])
    top.addBond(ace_atoms['C'],mapping[atom(residues[0],'N')])
    nhe=top.addResidue('NHE',chain,'15');na=top.addAtom('N',app.element.nitrogen,nhe);coords.append(amiden)
    top.addBond(mapping[atom(last,'C')],na)
    return app.Modeller(top,unit.Quantity(np.asarray(coords),unit.angstrom))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/2mi1_reduced_capped_ph75_20260924')
    parser.add_argument('--ph', type=float, default=7.5)
    parser.add_argument('--temperature', type=float, default=300.0)
    args = parser.parse_args()
    out = Path(args.output)
    if not out.is_absolute():
        out = ROOT/out
    out.mkdir(parents=True, exist_ok=False)
    pdb_path=ROOT/'data/structures/test.pdb'
    original_hash=core.sha(pdb_path)
    pdb=app.PDBFile(str(pdb_path))
    original_disulfides=[[a.residue.id,b.residue.id] for a,b in pdb.topology.bonds() if a.name==b.name=='SG']
    assert original_disulfides==[['3','14']],original_disulfides
    platform=mm.Platform.getPlatformByName('OpenCL')
    properties={'Precision':'mixed','DeviceIndex':'0','OpenCLPlatformIndex':'0'}
    protocol={'status':'running','openmm':mm.__version__,'pdbfixer':'installed 1.12.0; not used for capped preparation',
        'pdb_sha256':original_hash,'script_sha256':core.sha(__file__),
        'checkpoint_sha256':core.sha(ROOT/'models/site_predictor.pt'),
        'forcefield':['amber14-all.xml','amber14/tip3p.xml'],'platform':'OpenCL',
        'platform_properties':properties,'temperature_K':args.temperature,'preparation_pH':args.ph,
        'conditions_basis':'User confirmed N-acetyl, C-amide, unprotected Cys, TFA salt, pH 7.5 and 300 K. 0.15 M NaCl is an assumed solvent model; TFA/buffer concentrations unknown and TFA is not explicitly simulated.',
        'protonation_limit':'pH used for initial hydrogen assignment only; fixed protonation MD, not constant-pH MD. Both Cys explicitly modeled as neutral thiols; no site-specific pKa or thiolate population estimation.',
        'redox':'Remove SG-SG bond; neutral protonated CYS sidechains; classical fixed topology cannot reform disulfide',
        'input_inconsistency':'Original PDB has SG-SG connectivity and HG on both Cys. Discard all original hydrogens and explicitly rebuild reduced topology.',
        'original_disulfide_pairs':original_disulfides,
        'termini':'ACE (N-acetyl) and NHE (C-primary amide, NOT N-methylamide); cap heavy atoms geometrically initialized then minimized',
        'solvent':'TIP3P explicit water; 1 nm padding; NaCl 0.15 M plus neutralizing counterions',
        'ensemble':'NVT; no pressure coupling; fixed prepared solvent box',
        'nonbonded':'PME, 1 nm cutoff, default PME error tolerance',
        'constraints':'HBonds; default rigid water','timestep_fs':2,'friction_per_ps':1,
        'minimization_max_iterations':2000,'minimization_tolerance_kJ_mol_nm':10,
        'equilibration_ps':20,'production_ps_per_start':100,'snapshot_interval_ps':10,
        'starts':pdb.getNumFrames(),'seed_base':20260924,
        'sampling_limit':'Short exploratory relaxation; not converged reduced-state equilibrium ensemble; no enhanced sampling',
        'selection':'All 10 starts, all 10 fixed-time production snapshots per start; no score-based selection'}
    core.write_json(out/'protocol.json',protocol)
    forcefield=app.ForceField(*protocol['forcefield'])
    # NHE exists in ff14SB but is absent from the default hydrogen definitions.
    hydrogen_xml='<Residues><Residue name="NHE"><H name="HN1" parent="N"/><H name="HN2" parent="N"/></Residue></Residues>'
    app.Modeller.loadHydrogenDefinitions(io.StringIO(hydrogen_xml))
    snapshots, runs=[] ,[]
    started=time.perf_counter()
    for frame in range(pdb.getNumFrames()):
        index=frame+1; seed=20260924+index
        random.seed(seed);np.random.seed(seed)
        folder=out/f'start_{index:02d}';folder.mkdir()
        print(f'Start {index}/10: repair and reduction',flush=True)
        # Source standard-residue heavy atoms are complete. C-amidation needs no OXT.
        missing={};terminal={'ACE0':['CH3','C','O'],'NHE15':['N']}
        modeller,removed=no_disulfide(pdb.topology,pdb.getPositions(frame=frame))
        modeller=cap_termini(modeller)
        # PDBFixer may already omit the inconsistent source SG-SG connection.
        assert len(removed)<=1,removed
        assert not any(a.name==b.name=='SG' for a,b in modeller.topology.bonds())
        variants=['CYS' if r.name=='CYS' else None for r in modeller.topology.residues()]
        modeller.addHydrogens(forcefield,pH=args.ph,variants=variants,platform=mm.Platform.getPlatformByName('CPU'))
        peptide_topology=modeller.topology
        n_peptide=peptide_topology.getNumAtoms()
        for residue in peptide_topology.residues():
            if residue.name=='CYS':
                assert {'SG','HG'} <= {a.name for a in residue.atoms()}
            if residue.name=='NHE':
                assert {a.name for a in residue.atoms()}=={'N','HN1','HN2'}
        assert not any(a.name==b.name=='SG' for a,b in peptide_topology.bonds())
        write_pdb(folder/'reduced_prepared.pdb',peptide_topology,modeller.positions)
        ca_indices=[a.index for a in peptide_topology.atoms() if a.name=='CA']
        sg_indices=[a.index for a in peptide_topology.atoms() if a.name=='SG']
        reference=np.asarray(modeller.positions.value_in_unit(unit.angstrom))[ca_indices]
        modeller.addSolvent(forcefield,model='tip3p',padding=1*unit.nanometer,ionicStrength=0.15*unit.molar)
        system=forcefield.createSystem(modeller.topology,nonbondedMethod=app.PME,nonbondedCutoff=1*unit.nanometer,constraints=app.HBonds)
        integrator=mm.LangevinMiddleIntegrator(args.temperature*unit.kelvin,1/unit.picosecond,0.002*unit.picoseconds)
        integrator.setRandomNumberSeed(seed)
        simulation=app.Simulation(modeller.topology,system,integrator,platform,properties)
        simulation.context.setPositions(modeller.positions)
        (folder/'system.xml').write_text(mm.XmlSerializer.serialize(system),encoding='utf-8')
        (folder/'integrator.xml').write_text(mm.XmlSerializer.serialize(integrator),encoding='utf-8')
        write_pdb(folder/'solvated_initial.pdb',modeller.topology,modeller.positions)
        simulation.reporters.append(app.StateDataReporter(str(folder/'energy.csv'),1000,step=True,time=True,potentialEnergy=True,kineticEnergy=True,temperature=True))
        def snapshot(stage, time_ps):
            state=simulation.context.getState(getPositions=True,getEnergy=True,enforcePeriodicBox=False)
            energy=float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
            pos=state.getPositions(asNumpy=True)[:n_peptide]
            array=np.asarray(pos.value_in_unit(unit.angstrom))
            assert np.isfinite(array).all() and np.isfinite(energy)
            ca=array[ca_indices]; a=ca-ca.mean(axis=0); b=reference-reference.mean(axis=0)
            u,_,vh=np.linalg.svd(a.T@b); adjust=np.eye(3);adjust[2,2]=np.linalg.det(u@vh)
            rmsd=float(np.sqrt(np.mean(np.sum((a@(u@adjust@vh)-b)**2,axis=1))))
            path=folder/f'{stage}_{time_ps:03d}ps.pdb'
            write_pdb(path,peptide_topology,pos)
            row={'start':index,'stage':stage,'time_ps':time_ps,'pdb':path.relative_to(out).as_posix(),
                'sg_distance_A':float(np.linalg.norm(array[sg_indices[0]]-array[sg_indices[1]])),
                'ca_rmsd_from_prepared_A':rmsd,'solvated_potential_energy_kJ_mol':energy,'sha256':core.sha(path)}
            snapshots.append(row)
            core.write_jsonl(out/'snapshots.jsonl',snapshots)
        simulation.minimizeEnergy(tolerance=10*unit.kilojoule_per_mole/unit.nanometer,maxIterations=2000)
        snapshot('minimized',0)
        simulation.context.setVelocitiesToTemperature(args.temperature*unit.kelvin,seed)
        simulation.step(10000)
        snapshot('equilibrated',20)
        for step in range(1,11):
            simulation.step(5000)
            snapshot('production',20+step*10)
        simulation.saveState(str(folder/'final_state.xml'))
        simulation.saveCheckpoint(str(folder/'final_checkpoint.chk'))
        run={'start':index,'seed':seed,'missing_atoms_repaired':missing,'missing_terminal_atoms':terminal,
            'removed_disulfide_pairs':removed,'peptide_atoms':n_peptide,'solvated_atoms':modeller.topology.getNumAtoms(),
            'device':platform.getPropertyValue(simulation.context,'DeviceName'),'status':'complete'}
        runs.append(run);core.write_json(out/'runs.json',runs)
        print(f'Start {index}/10 complete; elapsed {time.perf_counter()-started:.1f}s; SG {snapshots[-1]["sg_distance_A"]:.2f} A',flush=True)
        del simulation,integrator,system
    assert core.sha(pdb_path)==original_hash
    assert core.sha(ROOT/'models/site_predictor.pt')==protocol['checkpoint_sha256']
    protocol.update(status='complete',elapsed_seconds=time.perf_counter()-started,production_snapshots=sum(r['stage']=='production' for r in snapshots))
    core.write_json(out/'protocol.json',protocol)
    print('MD COMPLETE',flush=True)


if __name__=='__main__':
    main()
