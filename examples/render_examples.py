"""Generate the example renders shown in the README.

Run from the repository root:

    python examples/render_examples.py

Writes an SVG and a PNG for each structure and reaction into this folder.
Everything here is drawn by the renderer itself.
"""
import os
import sys

# Make the top-level modules importable when run from the repo root.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import structure_draw as draw
import structure_svg as svg
import reaction_render as rxn

HERE = os.path.dirname(os.path.abspath(__file__))


def out(name):
    return os.path.join(HERE, name)


def structure(name, smiles):
    """Write <name>.svg and <name>.png for a single structure."""
    markup, _, _ = svg.render_svg(smiles)
    with open(out(f"{name}.svg"), "w") as f:
        f.write(markup)
    draw.render(smiles, out(f"{name}.png"))
    print("structure:", name)


def reaction(name, reactants, products, reagent="", conditions="", arrow="->"):
    markup, _, _ = rxn.render_reaction_svg(
        reactants, products, reagent=reagent, conditions=conditions, arrow=arrow
    )
    with open(out(f"{name}.svg"), "w") as f:
        f.write(markup)
    img = rxn.render_reaction_png(
        reactants, products, reagent=reagent, conditions=conditions, arrow=arrow
    )
    img.save(out(f"{name}.png"))
    print("reaction:  ", name)


if __name__ == "__main__":
    # Structures: a chain acid, an ester, and a substituted benzene.
    structure("ethanoic-acid", "CC(=O)O")
    structure("ethyl-ethanoate", "CCOC(C)=O")
    structure("nitrobenzene", "O=[N+]([O-])c1ccccc1")

    # Reaction 1: esterification, drawn with an equilibrium arrow.
    # "$H2O" is a literal token; the "$" marks non-SMILES species and is not
    # printed. Reagent and condition text are plain formula text (no "$").
    reaction(
        "esterification",
        ["CCO", "CC(=O)O"],
        ["CCOC(C)=O", "$H2O"],
        reagent="conc. H2SO4",
        conditions="heat",
        arrow="<=>",
    )

    # Reaction 2: oxidation of ethanol, with a literal oxidant token "$2[O]".
    reaction(
        "ethanol-oxidation",
        ["CCO", "$2[O]"],
        ["CC(=O)O", "$H2O"],
        reagent="K2Cr2O7 / H2SO4",
        conditions="reflux",
        arrow="->",
    )

    print("done")
