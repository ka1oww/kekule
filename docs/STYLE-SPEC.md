# A-level organic drawing: style specification

The single source of truth for how this renderer draws structures, and the rubric the evaluation
loop grades against. It exists so that anyone who installs the skill and asks for a structure gets
a drawing that reads as authentic A-level (H2 Chemistry 9729) house style, not a generic
cheminformatics depiction.

The engine is **deterministic rules over RDKit geometry**, not a learned model. So this spec is not
training data: it is the executable contract. Every rule below is either enforced in code, checked
by `tests/invariants.py`, or graded by a vision judge against real notes. Status tags:

- `[done]` implemented and invariant-checked
- `[partial]` implemented, edge cases open
- `[todo]` not yet handled

Evidence counts in the "Convention evidence" section are filled by the corpus-mining pass
(`corpus/`), which reads real A-level notes and question guides and records how each convention is
actually drawn, with page citations. Where the notes disagree with a rule here, the notes win and
this spec is amended.

---

## 1. The representations

| Form | Definition | Hydrogens | Carbons |
|------|-----------|-----------|---------|
| **Structural (condensed)** | groups written as text labels along a backbone; multiple bonds (C=O/C=C/C#C) drawn | absorbed into labels (`CH3`, `OH`, `NH2`) | shown as `C`/`CH`/`CH2`/`CH3` |
| **Natural (expanded)** | the **carbon/hetero BACKBONE drawn out on a 90-deg grid** (each backbone atom shown with its H at right angles), but every **substituent condensed to a label** (`CH3`, `OH`, `NH2`, `Br`, `OCH3`...); C=O / C=C drawn. Per the NJC naming-diagram style — "all bonds shown, real angles NOT followed, substituents condensed" | backbone H shown; substituent H in labels | backbone `C`; substituents condensed |
| **Displayed** | every atom and bond drawn at the **house geometry that DOES follow bond angle**: sp3 90-deg comb, sp2 120-deg trigonal, O/N bent | all shown, each on its own bond | shown as `C` |
| **Skeletal** | zig-zag lines; atoms implied at vertices/ends | implicit on C; shown only on heteroatoms | implicit at vertices |

Also: **stereochemical** (flying-wedge enantiomer pairs) and **geometric** (cis/trans C=C).

Structural is the NJC house default (what teaching material uses). **Displayed uses right-angle
"comb" geometry** by default (authentic A-level, per mined notes; `angles='natural'` gives the
tetrahedral variant). **Skeletal uses the ~120 degree zig-zag.** Rings in skeletal (and, for now,
displayed) render as the house hexagon + inscribed circle / plain polygon. `[done]`

## 2. Backbone and orientation
- Longest carbon chain laid **horizontal** (principal-axis rotation; natural angles preserved). `[done]`
- **Functional group oriented to the right**; a **left-hand terminal group reverses** so the bonded
  atom sits next to the bond and the rest hangs off the far side (`H3C-`, `H2N-`, `O2N-`, `HO-`).
  Driven by the `reversible` set + the renderer's bond-direction gate (`has_r and not has_l`), so it is
  now consistent across **all forms** — structural, natural, displayed and skeletal (skeletal/structural
  previously returned an empty set and left `-OH`/`-NH2` unflipped on the left end). `[done]`
- Molecule mirrored on x when needed so the functional-weighted centroid sits right of centre. `[done]`

## 3. Carbonyl group
- **C=O always points straight up** in structural; up-biased in displayed/skeletal via a y-flip. `[done]`
- Double bond drawn slightly short in structural (0.82x). Full length in displayed/skeletal. `[partial]` (consistency open)
- Aldehyde carbonyl **H drawn explicitly** even in condensed form. `[done]`
- Amide drawn out as `-C(=O)-NH2`, not a merged label. `[done]` (structural)

## 4. Hydrogen treatment
- Condensed: H absorbed into the heavy-atom label. `[done]`
- Displayed: every H is its own atom on a shortened bond (1.25 vs 1.95 heavy) so H's tuck in. `[done]`
- Skeletal: no H on carbon; H on O/N shown as part of the label (`OH`, `NH2`). An **in-chain
  secondary-amine N-H** (N with two heavy neighbours) is the exception: the N is labelled bare `N`
  and its H drawn stacked ABOVE it, so both backbone bonds meet the N, not the H. `[done]`

## 5. Bonds and geometry
- Single / double / triple as 1 / 2 / 3 parallel lines. `[done]`
- Bonds inset from a labelled atom so the line meets the letter, not its centre; blank vertices close
  fully (rings). `[done]`
- Stroke ~3px, bond length ~1.95 units (uniform: C-H = C-C = C-Br), monochrome black on white. `[done]`
- **Displayed hybridisation geometry** (unified BFS: every atom places its own neighbours off the bond
  it arrived on, by ITS hybridisation): `[done]`
  - **sp3 carbon** -> right-angle "comb": backbone bond continues straight, pendant H / branch bonds at 90 deg.
  - **sp2 carbon (any C=O, C=C, C=N)** -> trigonal ~120 deg splay; the =O is biased straight up.
  - **divalent O / N (O-H, N-H, ester/ether -O-, secondary amine C-N-C)** -> BENT ~120 deg, NEVER collinear.
  - **terminal cap (-CH3) on a non-horizontal bond** -> H's fan away from the parent (a hanging methyl
    would collide if combed onto the grid); horizontal chain-end methyls keep the clean comb.
  - Verified by an overlap/vision QA sweep across 30 acyclic molecules and by `tests/invariants.py`
    (`no_overlap`, `carbonyl_up`, `functional_right`, connectivity, bond count).

## 6. Aromatic and cyclic
- Benzene ring: **closed hexagon + inscribed circle** (delocalised), no C/H on the ring, in **all
  three forms**. `[done]`
- Cycloalkanes: plain polygons (no circle), upright; ring C=C as an inner parallel line. `[done]`
- **Ring substituents are form-aware** (`_ring_layout(mol, form)`): structural = condensed label;
  displayed = terminal `-CH3`/`-OH`/`-NH2` expanded to explicit atoms with the comb/bent geometry
  (O-H and N-H bend), and **carbonyl substituents (`-CHO`, `-COOH`, `-COCl`, ester) drawn OUT** as
  ring-C(=O)-X (=O biased up, trigonal splay, acid O-H bent) — like the acyclic displayed; skeletal =
  alkyl substituent drawn as bare zig-zag stub (first bond radial/straight, then zig-zag), heteroatom
  groups keep their label. `[done]` (longer plain-alkyl substituents stay condensed in displayed — edge, `[partial]`)
- Fused / polycyclic rings: out of scope; should fall back to RDKit depiction and flag. `[todo]`

## 7. Functional-group labels
- Canonical labels: `OH NH2 COOH CHO CN NO2 Cl Br`; ester `-COO-`; acyl chloride `-COCl`. `[partial]`
- Confirm from notes: `COOH` vs `CO2H`; `CHO` orientation; ester drawn as `-C(=O)-O-` vs `-COO-`.

## 8. Charges, ions, radicals
- Formal charge as a superscript token (`O^-` -> O with superscript minus). `[done]`
- Carbocation / carbanion / zwitterion (amino acid) layout. `[partial]`
- Free radical: the unpaired-electron dot. `[todo]`

## 9. Stereochemistry  (rebuilt from the mined Isomerism chapters; see `corpus/stereo_synthesis.json`)
- **Optical, single chiral centre**: FLYING-WEDGE (natta) template, tuned by overlaying the NJC figure.
  Letter `C` centre (no asterisk by default). Bond length `L=1.7`. House pose / angles: **the LARGEST
  group points straight UP** (90 deg, drawn 0.82 L, on the mirror axis); the next carbon is a plain bond
  DOWN-RIGHT (-42 deg); the hashed DASH (usually H) points DOWN-LEFT (197 deg, 0.85 L, into the plane);
  the bold WEDGE (heteroatom / functional group) points DOWN-LEFT (-122 deg, out of the plane). Carboxyl
  labelled `CO2H` (not `COOH`); a strongly-left simple group reverses (`CH3`->`H3C`, `OH`->`HO`). The
  "main group up" pose is the house convention we fixed by overlay — do NOT drop it to a small-group-up
  variant. Verified across lactic acid / alanine / butan-2-ol / 2-chloropropanoic acid. `[done]`
- **Enantiomer pair**: two flying-wedge structures split by a vertical dashed MIRROR-PLANE line
  (labelled), captioned "enantiomers (non-superimposable mirror images)"; the mirror is a single
  x-reflection. `render_stereo_pair`. `[done]`
- **cis/trans geometric (C=C)**: horizontal explicit `C=C`, four plain substituents ~120 deg,
  reference groups same-side (cis) / diagonal (trans), labelled `cis`/`trans`. `render_geometric_pair`.
  `[done]`
- **Labelling policy**: NEVER emit R/S, E/Z, D/L or (+)/(-) on a structure (out of 9729 scope). `[done]`
- Config-aware parity (draw a SPECIFIED @/@@ enantiomer), multiple stereocentres (2^n), meso
  (internal plane of symmetry), ring cis/trans (wedge/dash on a polygon), Fischer projection. `[todo]`

## 10. Beyond single molecules
- Delocalisation: dashed/partial bonds, partial charges (delta+/delta-). `[todo]`
- Polymer repeat unit: brackets with `n` and trailing bonds. `[todo]`
- Mechanism arrows: full (two-electron) and half (single-electron) curly arrows, lone pairs.
  `[partial]` (reaction_render draws schemes; curly arrows todo)

## 11. Typography
- Sans-serif (Arial), monochrome, subscript digits dropped low, superscript charges raised. `[done]`
- Portability: bundle or gracefully fall back on the font so non-macOS installs work. `[todo]`

---

## Convention evidence (mined from real notes)

Mined from 11 real A-level sources (intro chapter, Alkene/Arene/Alcohol notes, memory sheet, four
CHEM-IS-TRY question banks, two summary sheets): **270 catalogued structures**, conventions ranked
by independent-source support. Full data in `corpus/catalog.json` + `corpus/synthesis.json`; the
prioritised worklist is in `docs/REFINEMENT-BACKLOG.md`. Tally: **6 conventions the engine already
matches, 7 partial, 11 not yet handled.**

Confirmed by the notes (engine already correct): per-form hydrogen treatment, horizontal backbone
with functional group right, sans-serif monochrome typography with subscripts, double/triple bond
as parallel lines, amide drawn out, aldehyde H explicit.

Corrections the notes forced (notes win over the spec):
- **Displayed formula uses right-angle "comb" geometry** (horizontal C-C spine, vertical C-H bonds),
  NOT tetrahedral zig-zag. 8 sources. The zig-zag (~120 degrees) belongs to the **skeletal** form
  only. This reverses the earlier "natural angles" choice for displayed; see the amendment note.
- Benzene = hexagon + inscribed circle in **all three** forms (currently structural only).
- Cycloalkanes = plain polygons, no circle; ring C=C as an inner parallel line.
- Chlorine typeset with an **italic lowercase l** everywhere (house quirk vs the digit 1).
- Acid label context-selectable (COOH / CO2H / reversed HOOC- from the left); displayed acid = -C(=O)-OH.
- Amino acids are the backbone-orientation exception: NH2 left, COOH right; ionisation state matches conditions.
- Stereochemistry (wedge/dash, optical mirror pairs, cis/trans) is a formally defined A-level
  convention and completely unhandled.
