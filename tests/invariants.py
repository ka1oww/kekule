"""Property-based invariants for the layout engine.

No golden images, no labelled data: each check asserts a geometric/structural property that a
*correct* A-level drawing must satisfy, then runs it across a molecule set in all three forms.
A human reviews only the failures. This is the cheap, scalable half of the eval loop; the golden
set + vision judge (see docs/STYLE-SPEC.md) is the other half.

Run:  python3 tests/invariants.py           # prints a report, exits non-zero on any FAIL
"""
import os, sys, math, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rdkit import Chem
import structure_draw as J

# molecule set spanning the acyclic functional groups the new forms target
MOLECULES = [
    ("methane", "C"), ("ethane", "CC"), ("propane", "CCC"), ("ethene", "C=C"),
    ("chloroethane", "CCCl"), ("2-chlorobutane", "CCC(Cl)C"),
    ("ethanol", "CCO"), ("propan-2-ol", "CC(O)C"),
    ("ethanal", "CC=O"), ("propanone", "CC(C)=O"),
    ("ethanoic acid", "CC(=O)O"), ("propanoic acid", "CCC(=O)O"),
    ("methyl ethanoate", "COC(C)=O"), ("ethyl ethanoate", "CCOC(C)=O"),
    ("ethanamide", "CC(N)=O"), ("ethylamine", "CCN"),
    ("ethanenitrile", "CC#N"), ("2-aminopropanoic acid", "CC(N)C(=O)O"),
    ("benzene", "c1ccccc1"), ("cyclohexane", "C1CCCCC1"),
    ("phenol", "Oc1ccccc1"), ("methylbenzene", "Cc1ccccc1"),
]
FORMS = ("structural", "displayed", "skeletal")
MINSEP = 0.55          # min distance (coord units) between two labelled atom centres
CARBONYL_TOL = 0.15    # =O must sit at least this far ABOVE its carbon (units) to count as "up"


def _labelled(atoms):
    return {i: (x, y) for i, (lab, x, y) in atoms.items() if lab}


def chk_finite(atoms, bonds, mol, form):
    bad = [i for i, (lab, x, y) in atoms.items() if not (math.isfinite(x) and math.isfinite(y))]
    return (not bad, f"non-finite coords at {bad}" if bad else "")


def chk_placed(atoms, bonds, mol, form):
    missing = {i for (a, b, o) in bonds for i in (a, b) if i not in atoms}
    return (not missing, f"bond endpoints missing from atoms: {missing}" if missing else "")


def chk_connected(atoms, bonds, mol, form):
    if not atoms:
        return (True, "")
    adj = collections.defaultdict(set)
    for a, b, o in bonds:
        adj[a].add(b); adj[b].add(a)
    seen, stack = set(), [next(iter(atoms))]
    while stack:
        u = stack.pop()
        if u in seen:
            continue
        seen.add(u); stack.extend(adj[u] - seen)
    missing = set(atoms) - seen
    return (not missing, f"{len(missing)} atom(s) disconnected from the graph" if missing else "")


def chk_no_overlap(atoms, bonds, mol, form):
    L = list(_labelled(atoms).items())
    worst = None
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            d = math.hypot(L[i][1][0] - L[j][1][0], L[i][1][1] - L[j][1][1])
            if worst is None or d < worst[0]:
                worst = (d, L[i][0], L[j][0])
    if worst and worst[0] < MINSEP:
        return (False, f"labels {worst[1]} & {worst[2]} only {worst[0]:.2f}u apart (<{MINSEP})")
    return (True, "")


def chk_functional_right(atoms, bonds, mol, form):
    if form == "displayed":
        return (True, "")   # displayed draws hydroxyls/amines DOWN as substituents, not on the right
    # functional-weighted centroid should sit right of (or at) the geometric centroid
    heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and a.GetIdx() in atoms]
    if mol.GetRingInfo().NumRings() > 0 or len(heavy) < 2:
        return (True, "")  # rule is defined for acyclic backbones
    xs = {a.GetIdx(): atoms[a.GetIdx()][1] for a in heavy}
    xg = sum(xs.values()) / len(xs)
    w = {a.GetIdx(): J._fweight(a) for a in heavy}
    if not sum(w.values()):
        return (True, "")  # no functional group (e.g. alkane) -> nothing to orient
    xf = sum(w[i] * xs[i] for i in xs) / sum(w.values())
    return (xf >= xg - 0.3, f"functional centroid x={xf:.2f} left of geometric x={xg:.2f}")


def chk_carbonyl_up(atoms, bonds, mol, form):
    for o in J._carbonyl_Os(mol):
        if o not in atoms:
            continue
        c = mol.GetAtomWithIdx(o).GetNeighbors()[0].GetIdx()
        if c in atoms and atoms[o][2] < atoms[c][2] + CARBONYL_TOL:
            return (False, f"carbonyl O idx {o} not above its carbon")
    return (True, "")


def chk_bondcount(atoms, bonds, mol, form):
    if mol.GetRingInfo().NumRings() > 0:
        return (True, "")  # ring layouts condense substituents; a raw bond count is not comparable
    if form == "displayed":
        exp = Chem.AddHs(mol).GetNumBonds()
    elif form == "skeletal":
        exp = mol.GetNumBonds()
    else:
        return (True, "")  # condensed structural absorbs/adds H selectively; not a clean count
    return (len(bonds) == exp, f"{len(bonds)} bonds drawn, expected {exp}" if len(bonds) != exp else "")


CHECKS = [chk_finite, chk_placed, chk_connected, chk_no_overlap,
          chk_functional_right, chk_carbonyl_up, chk_bondcount]


def run():
    fails = []
    print(f"{'molecule':<24}{'form':<12}{'result'}")
    print("-" * 60)
    for name, smi in MOLECULES:
        for form in FORMS:
            try:
                a, b, c, rev = J.layout(smi, form)
                mol = Chem.MolFromSmiles(smi)
                probs = []
                for chk in CHECKS:
                    ok, msg = chk(a, b, mol, form)
                    if not ok:
                        probs.append(f"{chk.__name__[4:]}: {msg}")
                if probs:
                    fails.append((name, form, probs))
                    print(f"{name:<24}{form:<12}FAIL  {probs[0]}")
                    for p in probs[1:]:
                        print(f"{'':<36}{p}")
                else:
                    print(f"{name:<24}{form:<12}ok")
            except Exception as e:
                fails.append((name, form, [f"exception: {e}"]))
                print(f"{name:<24}{form:<12}ERROR {e}")
    print("-" * 60)
    total = len(MOLECULES) * len(FORMS)
    print(f"{total - len(fails)}/{total} passed, {len(fails)} failed")
    return len(fails)


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
