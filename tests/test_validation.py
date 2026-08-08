"""Regression tests for Kekule's public SMILES validation boundary."""

import os
import sys
import unittest
from unittest import mock

from PIL import Image

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import reaction_render as reaction
import embed_docx
import structure_draw as draw
import structure_svg as svg


class ValidationBoundaryTests(unittest.TestCase):
    def test_valid_core_input_returns_a_reusable_molecule(self):
        mol = draw.validate_input("CC(=O)O")
        self.assertEqual(mol.GetNumAtoms(), 4)
        self.assertEqual(mol.GetNumBonds(), 3)
        for form in ("structural", "displayed", "skeletal"):
            atoms, bonds, _, _ = draw._layout_mol(mol, form)
            self.assertTrue(atoms, form)
            self.assertTrue(bonds, form)

    def test_chemical_name_is_not_treated_as_smiles(self):
        with self.assertRaisesRegex(draw.InvalidSmilesError, "not chemical names"):
            draw.layout("ethanoic acid")

    def test_malformed_smiles_has_a_typed_human_readable_error(self):
        with self.assertRaisesRegex(draw.InvalidSmilesError, "Invalid SMILES"):
            draw.render("not-smiles")

    def test_empty_or_non_string_input_is_invalid(self):
        for value in ("", "   ", None, 42):
            with self.subTest(value=value):
                with self.assertRaises(draw.InvalidSmilesError):
                    draw.validate_input(value)

    def test_disconnected_structure_is_rejected(self):
        with self.assertRaisesRegex(draw.DisconnectedStructureError, "separate list entries"):
            draw.layout("CCO.O")

    def test_salt_is_rejected_as_one_structure(self):
        with self.assertRaisesRegex(draw.DisconnectedStructureError, "2 fragments"):
            draw.render("[Na+].[Cl-]")

    def test_radical_is_rejected_until_the_electron_dot_is_supported(self):
        with self.assertRaisesRegex(draw.UnsupportedRadicalError, "unpaired-electron dots"):
            svg.render_svg("[CH3]")

    def test_every_unsupported_multi_ring_topology_is_rejected(self):
        cases = {
            "naphthalene": "c1ccc2ccccc2c1",
            "biphenyl": "c1ccccc1-c2ccccc2",
            "decalin": "C1CCC2CCCCC2C1",
            "norbornane": "C1CC2CCC1C2",
            "spiro": "C1CCC2(CC1)CCCC2",
        }
        for name, smiles in cases.items():
            with self.subTest(name=name):
                with self.assertRaisesRegex(draw.UnsupportedTopologyError, "at most one simple ring"):
                    draw.layout(smiles)

    def test_multi_centre_stereo_pair_is_rejected(self):
        smiles = "O=C(O)[C@H](O)[C@H](O)C(=O)O"
        with self.assertRaisesRegex(draw.UnsupportedStereochemistryError, "exactly one"):
            draw.render_stereo_pair(smiles)

    def test_supported_stereo_and_geometric_pairs_still_render(self):
        stereo = draw.render_stereo_pair("CC(O)C(=O)O")
        geometric = draw.render_geometric_pair("CC=CC")
        images = [stereo] + [image for _, image in geometric]
        for image in images:
            with self.subTest(size=image.size):
                self.assertGreater(image.width, 0)
                self.assertGreater(image.height, 0)

    def test_non_stereogenic_pair_requests_have_typed_errors(self):
        with self.assertRaises(draw.UnsupportedStereochemistryError):
            draw.render_stereo_pair("CC")
        with self.assertRaisesRegex(draw.UnsupportedStereochemistryError, "two different groups"):
            draw.render_geometric_pair("C=C")

    def test_stereo_svg_has_a_typed_unsupported_error(self):
        for call in (svg.render_svg, svg.render_svg_inner):
            with self.subTest(call=call):
                with self.assertRaisesRegex(draw.UnsupportedStereochemistryError, "SVG backend"):
                    call("CC(O)C(=O)O", form="stereo")

    def test_supported_simple_rings_still_render(self):
        for smiles in ("c1ccccc1", "C1CCCCC1", "C1=CCCCC1"):
            for form in ("structural", "displayed", "skeletal"):
                with self.subTest(smiles=smiles, form=form):
                    markup, width, height = svg.render_svg(smiles, form=form)
                    self.assertTrue(markup.startswith("<svg"))
                    self.assertGreater(width, 0)
                    self.assertGreater(height, 0)


class PublicEntryPointTests(unittest.TestCase):
    def test_structure_entry_points_share_the_typed_boundary(self):
        entry_points = (
            lambda: draw.layout("ethanoic acid"),
            lambda: draw.render("ethanoic acid"),
            lambda: svg.render_svg("ethanoic acid"),
            lambda: svg.render_svg_inner("ethanoic acid"),
        )
        for call in entry_points:
            with self.subTest(call=call):
                with self.assertRaises(draw.InvalidSmilesError):
                    call()

    def test_reaction_smiles_species_use_the_same_boundary(self):
        entry_points = (
            lambda: reaction.render_reaction_svg(["CCO.O"], ["CC=O"]),
            lambda: reaction.render_reaction_png(["CCO.O"], ["CC=O"]),
        )
        for call in entry_points:
            with self.subTest(call=call):
                with self.assertRaisesRegex(draw.DisconnectedStructureError, "separate list entries"):
                    call()

    def test_reaction_rejects_unsupported_topology(self):
        with self.assertRaises(draw.UnsupportedTopologyError):
            reaction.render_reaction_svg(["c1ccc2ccccc2c1"], ["c1ccccc1"])

    def test_docx_structure_entry_point_uses_the_same_boundary(self):
        from docx import Document

        doc = Document()
        run = doc.add_paragraph().add_run()
        with self.assertRaises(draw.InvalidSmilesError):
            embed_docx.add_structure(doc, run, "ethanoic acid")

    def test_reaction_accepts_separate_species_and_explicit_literals(self):
        markup, width, height = reaction.render_reaction_svg(
            ["CCO", "$H2O"], ["CC=O", "$H2"], reagent="Cu", conditions="heat"
        )
        self.assertTrue(markup.startswith("<svg"))
        self.assertGreater(width, 0)
        self.assertGreater(height, 0)

    def test_structure_render_parses_smiles_once(self):
        parser = draw.Chem.MolFromSmiles
        with mock.patch.object(draw.Chem, "MolFromSmiles", wraps=parser) as wrapped:
            draw.render("CC(=O)O")
        self.assertEqual(wrapped.call_count, 1)

    def test_reaction_parses_each_structural_species_once(self):
        parser = draw.Chem.MolFromSmiles
        with mock.patch.object(draw.Chem, "MolFromSmiles", wraps=parser) as wrapped:
            reaction.render_reaction_svg(["CCO", "$2[O]"], ["CC=O", "$H2O"])
        self.assertEqual(wrapped.call_count, 2)


class SupportedOutputStabilityTests(unittest.TestCase):
    STRUCTURES = {
        "ethanoic-acid": "CC(=O)O",
        "ethyl-ethanoate": "CCOC(C)=O",
        "nitrobenzene": "O=[N+]([O-])c1ccccc1",
    }
    REACTIONS = {
        "esterification": {
            "reactants": ["CCO", "CC(=O)O"],
            "products": ["CCOC(C)=O", "$H2O"],
            "reagent": "conc. H2SO4",
            "conditions": "heat",
            "arrow": "<=>",
        },
        "ethanol-oxidation": {
            "reactants": ["CCO", "$2[O]"],
            "products": ["CC(=O)O", "$H2O"],
            "reagent": "K2Cr2O7 / H2SO4",
            "conditions": "reflux",
            "arrow": "->",
        },
    }

    def test_tracked_structure_examples_are_unchanged(self):
        for name, smiles in self.STRUCTURES.items():
            with self.subTest(name=name):
                markup, _, _ = svg.render_svg(smiles)
                with open(os.path.join(ROOT, "examples", f"{name}.svg")) as fh:
                    self.assertEqual(markup, fh.read())
                actual = draw.render(smiles)
                with Image.open(os.path.join(ROOT, "examples", f"{name}.png")) as tracked:
                    expected = tracked.convert("RGB")
                self.assertEqual(actual.size, expected.size)
                self.assertEqual(actual.tobytes(), expected.tobytes())

    def test_tracked_reaction_examples_are_unchanged(self):
        for name, kwargs in self.REACTIONS.items():
            with self.subTest(name=name):
                markup, _, _ = reaction.render_reaction_svg(**kwargs)
                with open(os.path.join(ROOT, "examples", f"{name}.svg")) as fh:
                    self.assertEqual(markup, fh.read())
                actual = reaction.render_reaction_png(**kwargs)
                with Image.open(os.path.join(ROOT, "examples", f"{name}.png")) as tracked:
                    expected = tracked.convert("RGB")
                self.assertEqual(actual.size, expected.size)
                self.assertEqual(actual.tobytes(), expected.tobytes())


if __name__ == "__main__":
    unittest.main()
