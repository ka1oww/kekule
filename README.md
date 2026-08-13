# SMILES-to-vector renderer

Turn a SMILES string into a sharp vector drawing of a chemical structure, or compose a whole reaction scheme, as SVG (with a PNG fallback). RDKit parses the SMILES; a custom layout pass draws it in the condensed structural style used in teaching and reference material, rather than the fully expanded skeletal form a generic renderer produces.

## Examples

Every image below was produced by the renderer itself. Regenerate them with `python3 examples/render_examples.py`.

**Structures.** Chains, common functional groups, and substituted benzene rings:

![ethanoic acid](examples/ethanoic-acid.png)
![ethyl ethanoate](examples/ethyl-ethanoate.png)
![nitrobenzene](examples/nitrobenzene.png)

**Reactions.** Multiple species, reagent above and conditions below the arrow, and an equilibrium arrow where the chemistry calls for one:

![esterification](examples/esterification.png)

![oxidation of ethanol](examples/ethanol-oxidation.png)

## What it does

- **SMILES to vector.** Parses with RDKit, then lays the molecule out with a horizontal backbone, carbonyls drawn upward, and every hydrogen shown, matching the condensed convention.
- **Sharp SVG output.** Vector, so it stays crisp at any size. A raster PNG fallback is available through the same call.
- **Reaction schemes.** Composes several species into one drawing with reagent and condition text set over the arrow, straight or equilibrium arrows, formal charges, and built-up fractions for half-equation coefficients.
- **Direct Word export.** Writes the vector drawing straight into a `.docx`. See the note below on how, because it is the part I am most pleased with.

## Install

```bash
pip install -r requirements.txt
```

Requirements are `rdkit` and `pillow`, plus `python-docx` for Word export. CI supports Python 3.11 through 3.13. The renderer uses Arial when it is installed and falls back to Liberation Sans or DejaVu Sans on other platforms. Reference-image byte comparisons run only with Arial; the cross-platform gate asserts stable chemical and artifact semantics instead of font-dependent pixels.

## Usage

```python
import structure_svg as svg
markup, width, height = svg.render_svg("CC(=O)O")     # ethanoic acid, returns an SVG string

import structure_draw as draw
draw.render("CC(=O)O", "ethanoic-acid.png")           # write a PNG instead
draw.render("[NH3]->[BH3]", "ammonia-borane.png", dative=True)   # coordinate bond as an arrow

import reaction_render as rxn
markup, width, height = rxn.render_reaction_svg(
    ["CCO", "CC(=O)O"], ["CCOC(C)=O", "$H2O"],
    reagent="conc. H2SO4", conditions="heat", arrow="<=>",
)

from embed_docx import build_doc
build_doc([("ethanoic acid", "CC(=O)O"), ("benzene", "c1ccccc1")], "structures.docx")
```

In a reaction, most species are SMILES. Ionic compounds may use dot-disconnected SMILES such as `CC(=O)[O-].[Na+]`. Single-atom fragment tokens such as `[Na+]`, `[Cl-]`, and `[NH4+]` retain their labels and formal charges. A species that is not structural SMILES, an inorganic reagent, a radical, or a coefficient, is written as a literal token with a leading `$`, for example `$H2O`, `$2[O]`, `$Na^+`, `$conc. HCl`. Digits after a letter auto-subscript, so `$H2SO4` prints as H₂SO₄. The `$` marks the token and is not drawn. Reagent and condition text are plain formula text and take no `$`.

## Input validation

Every public structure and reaction entry point validates a structural species before layout. Kekule accepts simple acyclic or single-ring SMILES. A dot-disconnected request is accepted when every fragment meets that scope and has a visible atom label, bond, or circle in the selected representation. It does not convert chemical names to SMILES.

Invalid or unsupported input raises a typed `ValueError` subclass from `structure_draw`:

- `InvalidSmilesError` for malformed SMILES or a chemical name passed as SMILES.
- `DisconnectedStructureError` when a specialised renderer that needs one connected fragment, such as `stereo` or `geometric`, receives dot-disconnected SMILES, or when any fragment has no visible atom label, bond, or circle in the selected representation.
- `UnsupportedTopologyError` for multi-ring, fused, bridged, or spiro structures.
- `UnsupportedRadicalError` while structural radical-dot rendering is unsupported. Radicals may still be typeset explicitly as reaction literals such as `$•CH3`.
- `UnsupportedStereochemistryError` when a specialised stereo request cannot be represented faithfully. The raster stereo-pair helper supports exactly one tetrahedral stereocentre and draws both enantiomers; it does not select one configuration from `@` or `@@`. Multi-centre stereo is rejected, and wedge/dash stereo SVG is outside the supported scope.
- `UnsupportedDativeBondError` when SMILES dative notation is used without the `dative=True` opt-in, or when a dative arrow is requested outside the structural, displayed, and skeletal forms, or together with a ring. The message identifies the offending bond as donor and acceptor.

All of these inherit from `KekuleInputError`, so callers may catch either the specific condition or the shared base class.

## Scope

The layout engine targets acyclic molecules (chains with one functional group and simple branches) and simple single-ring structures, including single-benzene-ring derivatives. Dot-disconnected compounds are rendered fragment by fragment from left to right with a clear gap in displayed, structural, and skeletal forms when every fragment has a visible primitive. Each fragment retains the same scope checks, so fused, bridged, spiro, radical, and other unsupported fragments are rejected before layout and identified in the error.

A coordinate (dative) bond written with SMILES dative notation, such as the adduct `[NH3]->[BH3]`, is drawn as a single-headed arrow from the donor to the acceptor in the structural, displayed, and skeletal forms when the caller opts in with `dative=True` on `layout`, `render`, `render_svg`, or `render_svg_inner`. The default refuses dative input rather than drawing a plain line, so existing output is unchanged unless the option is requested. Formal charges are never reinterpreted as coordinate bonds: `[NH4+]` keeps its charged label, and an arrow appears only where the input actually expresses a dative bond.

## Readiness gate

Run the complete local preview gate with one command from the repository root:

```bash
python3 scripts/readiness_gate.py
```

This compiles the Python sources, runs the geometric invariant and input-validation suites, and executes every case in [`evaluation/preview_manifest.json`](evaluation/preview_manifest.json), including the portable DOCX package smoke. The manifest is the single human-readable inventory of supported preview examples and typed rejection cases; the runner contains only generic entry-point dispatch and semantic artifact checks. GitHub Actions installs the declared dependencies and invokes this same command on Python 3.11, 3.12, and 3.13.

## Licence

[MIT](LICENSE)
