# Refinement backlog (mined from real A-level notes)

Generated from 11 sources, 270 catalogued structures. Conventions ranked by number of independent sources that drew them that way. `MATCH` = engine already does this; `MISS`/`PARTIAL` = work to do. Full data in `corpus/`.

## Conventions by support

- **[10] MATCH** `hydrogen-treatment` — Hydrogen is shown per form: DISPLAYED draws every C-H on its own bond; SKELETAL shows no C-H (any undrawn 4th bond on a carbon is an implied C-H); CONDENSED folds H into CH3/CH2 group labels. Hydrogen on heteroatoms (O-H, N-H, NH2, NH3+) is ALWAYS drawn, even in skeletal.
  - -> Keep. Verify skeletal always renders heteroatom H (OH/NH/NH2) and add an optional 'defining C-H' flag (aldehyde H, amino-acid alpha C-H) which several sources draw explicitly even in skeletal.
- **[10] PARTIAL** `fg-labels` — Carboxylic acid label is inconsistent across and within sources: COOH, CO2H, and reversed HOOC-/HO2C- when the bond emerges from the left; the drawn/displayed form is -C(=O)-OH. A renderer must support all spellings, context-selectable.
  - -> Add an acid-label style option (COOH default, CO2H alt); auto-reverse to HOOC-/HO2C- when the group attaches from the left; render displayed acids as -C(=O)-OH.
- **[9] MATCH** `typography` — Atom labels sans-serif; digit counts as subscripts, charges/hybridisation as superscripts; structures monochrome black on white (colour only pedagogical: blue=answers, red=highlight/reagents/products). One source (RP Carbonyl) dissents with a serif face.
  - -> Keep Arial/sans-serif monochrome + subscripts. Optionally expose a serif toggle to match the Carbonyl RP source.
- **[8] MATCH** `backbone-orientation` — Carbon chain drawn horizontally left-to-right with the principal functional group at the RIGHT-hand end. Exception: amino acids are drawn NH2 on the left, COOH on the right.
  - -> Keep. Add the amino-acid exception (NH2 left / COOH right) and the memory-sheet rule that the amino-acid ionisation state matches reaction conditions (+H3N-...-COOH acidic vs H2N-...-COO- basic).
- **[8] MISS** `bond-geometry` — DISPLAYED formulae use RIGHT-ANGLE (90 degrees) 'comb' geometry: horizontal C-C spine with strictly vertical C-H/substituent bonds, NOT a tetrahedral zig-zag. SKELETAL formulae use ~120 degrees tetrahedral zig-zag; sp2 C=C substituents splay ~120 degrees.
  - -> HIGH: change displayed-form layout from 'natural tetrahedral/trigonal angles' to the orthogonal comb (horizontal spine, vertical H/substituent bonds). Keep 120-degree zig-zag for skeletal only.
- **[8] PARTIAL** `benzene` — Benzene and all derivatives are drawn as a regular hexagon with an INSCRIBED CIRCLE (delocalised) as the house default; ring H implicit.
  - -> HIGH: apply hexagon+inscribed-circle in ALL three forms; currently only the structural/condensed form is house-styled while displayed and skeletal rings are unstyled.
- **[8] MISS** `cycloalkanes` — Cycloalkanes/rings drawn as plain regular polygons with NO inscribed circle (this is how cyclohexane is distinguished from benzene); ring unsaturation added as a second parallel line on the affected edge (cyclohexene = 1, cyclohexadiene = 2). Substituents attach at a vertex.
  - -> HIGH: implement ring rendering (regular polygons, vertex attachment, inner second line for ring C=C) and honour the plain-polygon-vs-circle distinction.
- **[8] MATCH** `double-bond-drawing` — Double bond = two parallel lines; triple bond = three parallel lines; roughly full length on chains (some sources draw the inner line shorter/inset for ring C=C), no large deliberate shortening.
  - -> Keep. Optionally draw the second line shorter/inset for ring double bonds to match RP style.
- **[8] PARTIAL** `carbonyl-direction` — C=O UP is the modal default (notes almost always put O above the carbonyl C), but orientation is PRAGMATIC: terminal aldehydes often draw =O to the side, esters sometimes down, ring ketones up, and the two -COOH of a symmetric diacid point in OPPOSITE directions; one summary source draws -COOH sideways with O-H hanging down.
  - -> Keep C=O-up as default but allow per-group orientation overrides: terminal aldehyde to the side, ester down, opposite directions for diacids.
- **[8] MISS** `charges-ions` — Formal charges as superscript +/- next to the atom; carbanion/carboxylate sometimes a CIRCLED minus (with lone-pair dots); partial charges as delta+/delta-; protonated amine as left-superscript +H3N-; complex ions in square brackets with an outside superscript charge; lone pairs as dot pairs.
  - -> Add charge glyphs (superscript +/-, circled-minus option, delta+/-, +H3N- left-superscript, [ ]^n complex ions) and lone-pair dot pairs.
- **[7] MISS** `chlorine-typography` — Chlorine typeset 'Cl' with an ITALIC lowercase l everywhere (Cl2, HCl, AlCl3, CCl4, -Cl), a deliberate house quirk to distinguish l from digit 1; iodine uses a capital I; Br upright.
  - -> Cheap high-value fix: render the l in every chlorine label italic; keep I/Br upright.
- **[7] MISS** `reaction-scheme-arrows` — Reaction schemes: single straight arrow with reagent ABOVE and conditions BELOW; [O] for oxidation; equilibrium = double harpoon; reversible catalytic/aromatisation process = open double-headed arrow; resonance = double-headed straight arrow; crossed arrow / X = no reaction; memory sheet uses colour-coded paired forward/back arrows.
  - -> Separate scheme layer (out of single-structure scope): arrow primitives with above/below annotations, [O], and equilibrium/resonance/reversible/no-reaction variants.
- **[5] PARTIAL** `fg-depiction-ester` — Ester drawn with an EXPLICIT carbonyl -C(=O)-O- (never condensed to -COO- in the drawing); contracted only inline as -CO2-/-OCO-.
  - -> Render esters as -C(=O)-O-R with explicit carbonyl in all drawn forms; reserve -COO-/-CO2- for inline condensed text only.
- **[5] MISS** `fused-rings` — Fused ring systems share an edge and may mix polygon sizes (hexagon+pentagon indane/indene, two fused hexagons naphthalene/decalin, hexagon+square benzocyclobutene, benzo fused to N-ring); 3-membered epoxide = triangle with O at apex.
  - -> Add fused-ring support: shared-edge polygon fusion, mixed ring sizes, heteroatom ring members, epoxide triangle.
- **[5] PARTIAL** `label-ordering` — Condensed group labels are mirror-ordered by side: a group on the LEFT is written with its attaching atom last (H3C-, HO-, H2N-, O2N-, HOOC-); on the RIGHT with the attaching atom first (-CH3, -OH, -NH2, -NO2, -COOH), so the bonded atom always faces the bond.
  - -> Auto-reverse condensed group labels based on bond direction so the attaching atom touches the bond.
- **[4] MISS** `stereochemistry-wedge-dash` — Bold WEDGE = bond out of the plane (toward viewer); hashed/DASH = into the plane. Used for optical isomers as mirror-image PAIRS about a chiral C with four different groups, and for 3-D tetrahedral depictions. Definition stated verbatim in 10 Intro to Organic.
  - -> HIGH: add wedge/dash bond rendering and optical-isomer mirror-pair layout; a formally defined A-level convention that is completely unhandled.
- **[4] MATCH** `fg-depiction-amide` — Amide drawn -C(=O)-N(H)- with the N-H shown explicitly (also peptide bonds and anilides); not condensed to CONH2 in drawn form.
  - -> Keep (engine draws -C(=O)-NH2); ensure N-H is shown and support N-substituted (secondary) amides and peptide -C(=O)-NH-.
- **[4] PARTIAL** `aromatic-substituent-placement` — On the ring the principal/first-named group sits at the TOP vertex; ortho (1,2) = adjacent vertex, para (1,4) = opposite/bottom vertex, 2,4,6-trisubstitution fills the two upper-side vertices + bottom. (10 Intro uses a flat-side hexagon variant with para drawn horizontally left/right.)
  - -> Implement deterministic vertex placement (top / ortho-adjacent / para-opposite / 2,4,6) with pointy-top default; expose ring orientation.
- **[4] MISS** `mechanism-arrows` — Curly arrows flow from electron-rich to electron-poor; FULL-headed arrow = electron pair (heterolysis), HALF/fishhook = single electron (homolysis/radical); lone pairs as dot pairs; arenium/Wheland intermediate = ring with broken inscribed circle + '+' (or circled +).
  - -> Mechanism layer (separate feature): full/fishhook curly-arrow primitives, lone-pair dots, arenium broken-circle+plus rendering.
- **[4] MATCH** `fg-depiction-aldehyde` — Notes draw the aldehyde as a displayed -C(=O)-H cluster with H explicit (never text 'CHO' in the drawing); the RP arene bank instead condenses it to a 'CHO' label at the vertex. Genuine notes-vs-questionbank style conflict.
  - -> Keep explicit-H default (matches notes majority); add a 'CHO label' toggle for arene/question-bank style.
- **[3] MISS** `stereochemistry-cis-trans` — Geometric isomerism shown by placing the two substituents on the SAME side (cis/Z) or OPPOSITE sides (trans/E) of a horizontal C=C; several sources draw geometry but never label E/Z; many notes omit it entirely.
  - -> Add cis/trans substituent placement across C=C (same-/opposite-side) for geometric isomers.
- **[3] PARTIAL** `benzene-kekule` — Kekule ring (three alternating C=C, no circle) is used selectively: historical cyclohexa-1,3,5-triene, resonance canonical forms, NON-aromatic ring intermediates/products (arenium ions, Birch dienes, arene oxides), and as house style by some schools; the delocalised circle is the default everywhere else. Compilations can show both, even within one question.
  - -> Add a Kekule ring option for non-aromatic ring intermediates/resonance and school-specific styles; default stays circle.
- **[1] MISS** `free-radicals` — Free radical shown as a bold DOT on the atom (homolysis products of e.g. free-radical substitution).
  - -> Low: add radical-dot glyph (needed only for free-radical substitution).
- **[1] MISS** `polymer-repeat-unit` — Polymer repeat unit enclosed in square brackets with the two continuation bonds crossing the brackets and a subscript italic 'n'; an 'n' prefix on the monomer denotes polymerisation.
  - -> Low: add repeat-unit brackets + crossing continuation bonds + subscript n.

## Engine change backlog

### HIGH
- **House-style aromatic rings (hexagon + inscribed circle) in the displayed and skeletal forms, not only the condensed/structural form.** — Benzene = hexagon + inscribed circle is the single most-supported drawing convention (8 sources) yet the engine leaves rings unstyled in 2 of its 3 forms.
- **Switch displayed-form bond geometry to the right-angle orthogonal 'comb' (horizontal spine + vertical H/substituent bonds); keep the 120-degree zig-zag only for skeletal.** — 8 sources draw displayed formulae with right angles; the engine currently uses natural tetrahedral/trigonal angles, a direct mismatch with the house 'displayed' look.
- **Ring engine: render cycloalkanes/rings as plain regular polygons (no circle), attach substituents at vertices, draw ring C=C as a second inner edge line, pointy-top/vertex-up default.** — 8 sources; the plain-polygon-vs-circle distinction is how students tell cyclohexane from benzene, and rings are currently unstyled in displayed and skeletal.
- **Wedge/dash stereochemistry plus optical-isomer mirror-pair layout.** — A formally defined A-level convention (10 Intro) drawn in 4 sources; completely unhandled by the engine.

### MED
- **cis/trans geometry across C=C via same-side/opposite-side substituent placement.** — 3 sources draw geometric isomers; examinable and unhandled.
- **Italic lowercase l in every chlorine label (Cl); keep I/Br upright.** — Deliberate house quirk in 7 sources; a cheap typographic change with high recognisability value.
- **Carboxylic-acid label options (COOH default, CO2H alt) with automatic HOOC-/HO2C- reversal when attaching from the left; render displayed acids as -C(=O)-OH.** — Acid spelling varies in all 10 data-bearing sources; a deterministic renderer must pick and reverse correctly.
- **Context-aware carbonyl orientation: keep C=O-up default but allow terminal aldehyde =O to the side, ester down, and opposite directions for the two -COOH of a diacid.** — 8 sources show C=O direction is pragmatic; the engine hard-codes up.
- **Condensed-label mirror-ordering by bond side (H3C-/-CH3, HO-/-OH, H2N-/-NH2, O2N-/-NO2).** — 5 sources; the bonded atom must face the bond for correct, professional-looking labels.
- **Fused/polycyclic ring support: shared-edge fusion, mixed polygon sizes, heteroatom ring members, epoxide triangle.** — 5 sources (naphthalene, indane/indene, decalin, benzocyclobutene, tetrahydroquinoline, epoxides); unhandled.
- **Charge/ion glyphs: superscript +/-, circled-minus option for carbanion/carboxylate, delta+/- partial charges, left-superscript +H3N-, [ ]^n complex ions, lone-pair dot pairs.** — 8 sources; needed for anions/cations, reactive intermediates, and condition-dependent acid-base forms.
- **Deterministic aromatic substituent placement (principal group top vertex, ortho = adjacent, para = opposite/bottom, 2,4,6 pattern).** — 4 sources; substituted-arene rendering needs fixed, reproducible positions.
- **Draw esters with an explicit -C(=O)-O- carbonyl (not -COO-) in all drawn forms.** — 5 sources; contraction to -COO-/-CO2- is reserved for inline text only.

### LOW
- **Kekule ring option for non-aromatic ring intermediates/resonance forms and school-specific house styles.** — 3 sources use Kekule selectively (arenium ions, Birch dienes, arene oxides, resonance canonical forms).
- **Aldehyde CHO-label toggle (arene/question-bank style) alongside the explicit -C(=O)-H default.** — Notes draw -C(=O)-H (3 sources) but the RP arene bank labels it CHO, a genuine documented style conflict.
- **Reaction-scheme layer: single/equilibrium/reversible/resonance/no-reaction arrows with reagent-above and conditions-below, plus [O].** — 7 sources render whole schemes; out of single-structure scope but needed for full-equation output.
- **Mechanism primitives: full-headed and fishhook curly arrows, lone pairs, arenium broken-circle+plus.** — 4 sources; a large separate feature for mechanism drawing.
- **Free-radical bold-dot glyph.** — 1 source (free-radical substitution); niche but currently impossible.
- **Polymer repeat-unit brackets with crossing continuation bonds + subscript n.** — 1 source (polystyrene); niche addition-polymer need.
- **Optional serif label toggle.** — One source (RP Carbonyl) uses a serif face against the sans-serif majority; low-impact parity option.


## Batch audit vs source (100 compounds)

Rendered all forms for 100 mined compounds and diffed each against the real notes (grouped by source
PDF). Result after the atom-preserving condenser fix: **50 match / 40 close / 8 wrong / 2 not-found, avg 78** (up from avg 69 with the
old lossy condenser). Full data in `corpus/batch_audit.json`.

Remaining, by category:
- **Engine scope (single-ring only)**: fused/bicyclic rings (lactones, indanone, decalin) open into
  chains. Biggest build; out of the original design.
- **Real fixable**: displayed form of ring compounds reuses the condensed structural (not drawn out);
  aldehyde on a left ring-vertex reads `HOC` not `OHC`; charges drift in the skeletal pane.
- **NOT the renderer (mined-catalog data)**: cis/trans lost because SMILES lack `/\` marks; some
  catalog SMILES/names disagree (e.g. "3-chlorocyclohexanol" SMILES is the 1,4 isomer). Ring
  substituent *placement* itself is correct (1,2/1,3/1,4 verified distinct).
- **Cosmetic house-convention**: pointy-top rings, italic `l` in Cl, `CO2H` spelling.

## R2 diagnosis and rerun

The eight `WRONG` records were diagnosed before rule changes and fall into four classes. The
classification distinguishes form selection from geometry, and does not treat an unsupported
topology as a reason to widen the renderer.

| Class | Examples | Form-selection or geometry | Decision |
|---|---|---|---|
| Displayed form lost directional alkene geometry | trans-1,2-dichloroethene; trans-but-2-ene | Form-selection plus geometry. The displayed path ignored SMILES E/Z and made both isomers look cis. | Fixed by seating the C=C horizontally and placing substituents on source-backed trigonal sides. |
| Source record encoded a different compound | 3-chlorocyclohexanol; 1-methylcyclopenta-1,3-diene | Neither. The catalog SMILES disagreed with the named/source structure. | Fixed the catalog and audit records to `OC1CC(Cl)CCC1` and `CC1=CC=CC1`. |
| Ring-attached branch was collapsed in the selected form | Vanillin | Form-selection. Methoxy, aldehyde, and the mixed displayed branch were not expanded together. | Fixed displayed ring branches while retaining the one-ring scope. |
| Fused or bicyclic topology was linearised | Iodolactone, indan-1-one, decahydronaphthalen-4a-ol | Neither. These are outside the accepted single-ring scope. | Kept typed `UnsupportedTopologyError`; no fallback drawing and no scope expansion. |

Representative `CLOSE` records showed the same form-selection class in paracetamol, propyl
benzoate, benzaldehyde, styrene, and benzyl alcohol. Other close records were geometry or
typography details: skeletal carbocation charges, cis geometry, chlorine typography, and left-facing
acid/aldehyde label ordering. The rule fixes cover those classes where the backlog has evidence
support. Resonance pairs, source-specific acid spelling, ring orientation, missing contextual
annotations, and mixed-form source disagreements remain findings rather than forced conventions.

The repository does not contain the source-PDF crops or a visual adjudicator, so the original
100-compound `MATCH` / `CLOSE` / `WRONG` / `NOT FOUND` score cannot be honestly recomputed here.
The stored baseline remains **50 / 40 / 8 / 2, average 77.6**. The R2 rule-level replay is recorded
in `corpus/batch_audit.json`: it names the four formerly wrong records fixed by rules or corrected
source data, the one formerly wrong record fixed by form selection, and the three formerly wrong
records now rejected at the typed scope boundary. It reports no invented visual score. The
automated rerun passed **66/66 invariants, 27/27 validation tests, and 21/21 preview cases**.

No regressions were observed in the tracked example artifacts or in the typed rejection path. The
remaining visual score is an explicit follow-up for the source-backed adjudication environment,
not evidence to be replaced with a model or training step.
