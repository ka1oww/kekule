"""Execute the readable narrow-preview manifest as a semantic CI gate."""

import json
import math
import os
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile

from PIL import Image, ImageChops
from rdkit import Chem

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "tests"))

import embed_docx
import invariants
import reaction_render
import structure_draw
import structure_svg


MANIFEST_PATH = os.path.join(ROOT, "evaluation", "preview_manifest.json")


def load_manifest():
    with open(MANIFEST_PATH, encoding="utf-8") as handle:
        manifest = json.load(handle)
    ids = [case["id"] for case in manifest["cases"]]
    if len(ids) != len(set(ids)):
        raise ValueError("Preview manifest case identifiers must be unique")
    return manifest


def _call(case, temp_dir):
    entry_point = case["entry_point"]
    value = case["input"]
    representation = case["representation"]

    if entry_point == "structure_draw.validate_input":
        return structure_draw.validate_input(value, representation)
    if entry_point == "structure_draw.layout":
        return structure_draw.layout(value, representation)
    if entry_point == "structure_draw.render":
        return structure_draw.render(value, form=representation)
    if entry_point == "structure_draw.render_stereo_pair":
        return structure_draw.render_stereo_pair(value)
    if entry_point == "structure_draw.render_geometric_pair":
        return structure_draw.render_geometric_pair(value)
    if entry_point == "structure_svg.render_svg":
        return structure_svg.render_svg(value, form=representation)
    if entry_point == "reaction_render.render_reaction_svg":
        return reaction_render.render_reaction_svg(**value)
    if entry_point == "embed_docx.build_doc":
        output = os.path.join(temp_dir, f"{case['id']}.docx")
        embed_docx.build_doc(value, output, title="Kekule preview smoke")
        return output
    raise ValueError(f"Unknown manifest entry point: {entry_point}")


def _assert_source_molecule(test, case):
    expected = case["expected"].get("molecule")
    if expected is None:
        return None
    molecule = Chem.MolFromSmiles(case["input"])
    test.assertIsNotNone(molecule)
    actual = {
        "atoms": molecule.GetNumAtoms(),
        "bonds": molecule.GetNumBonds(),
        "rings": molecule.GetRingInfo().NumRings(),
        "formal_charge": Chem.GetFormalCharge(molecule),
    }
    test.assertEqual(actual, expected)
    return molecule


def _assert_nonblank_image(test, image):
    test.assertIsInstance(image, Image.Image)
    test.assertGreater(image.width, 0)
    test.assertGreater(image.height, 0)
    white = Image.new(image.mode, image.size, "white")
    test.assertIsNotNone(ImageChops.difference(image, white).getbbox())


def _assert_svg(test, result):
    markup, width, height = result
    root = ET.fromstring(markup)
    test.assertTrue(root.tag.endswith("svg"))
    test.assertGreater(width, 0)
    test.assertGreater(height, 0)
    test.assertEqual(root.attrib["width"], str(width))
    test.assertEqual(root.attrib["height"], str(height))


def _assert_reaction_svg(test, case, result):
    _assert_svg(test, result)
    markup = result[0]
    root = ET.fromstring(markup)
    groups = [child for child in root if child.tag.endswith("g")]
    reactants = case["input"]["reactants"]
    products = case["input"]["products"]
    expected_blocks = 2 * (len(reactants) + len(products)) - 1
    test.assertEqual(len(groups), expected_blocks)
    arrow_group = groups[2 * len(reactants) - 1]
    arrow_elements = [child.tag.rsplit("}", 1)[-1] for child in arrow_group]
    if case["input"].get("arrow", "->") == "<=>":
        test.assertEqual(arrow_elements.count("line"), 4)
        test.assertNotIn("path", arrow_elements)
    else:
        test.assertEqual(arrow_elements.count("line"), 1)
        test.assertEqual(arrow_elements.count("path"), 1)


def _assert_connected_layout(test, result):
    atoms, bonds = result[0], result[1]
    test.assertTrue(atoms)
    test.assertTrue(bonds)
    for _, x, y in atoms.values():
        test.assertTrue(math.isfinite(x) and math.isfinite(y))
    adjacency = {atom_id: set() for atom_id in atoms}
    for first, second, _ in bonds:
        test.assertIn(first, atoms)
        test.assertIn(second, atoms)
        adjacency[first].add(second)
        adjacency[second].add(first)
    seen = set()
    pending = [next(iter(atoms))]
    while pending:
        atom_id = pending.pop()
        if atom_id not in seen:
            seen.add(atom_id)
            pending.extend(adjacency[atom_id] - seen)
    test.assertEqual(seen, set(atoms))


def _assert_docx_package(test, path, structure_count):
    test.assertTrue(zipfile.is_zipfile(path))
    with zipfile.ZipFile(path) as package:
        names = set(package.namelist())
        test.assertIn("[Content_Types].xml", names)
        test.assertIn("word/document.xml", names)
        test.assertEqual(sum(name.endswith(".png") for name in names), structure_count)
        test.assertEqual(sum(name.endswith(".svg") for name in names), structure_count)
        document_xml = package.read("word/document.xml")
        test.assertEqual(document_xml.count(b"svgBlip"), structure_count)


def _assert_structure_invariants(test, case, molecule, result):
    representation = case["representation"]
    if molecule is None or representation not in invariants.FORMS:
        return
    if case["expected"]["artifact"] == "connected_layout":
        layout = result
    else:
        layout = structure_draw.layout(case["input"], representation)
    atoms, bonds = layout[0], layout[1]
    for check in invariants.CHECKS:
        ok, message = check(atoms, bonds, molecule, representation)
        test.assertTrue(ok, f"{check.__name__}: {message}")


def _assert_supported(test, case, result):
    molecule = _assert_source_molecule(test, case)
    _assert_structure_invariants(test, case, molecule, result)
    artifact = case["expected"]["artifact"]
    if artifact == "svg":
        _assert_svg(test, result)
    elif artifact == "reaction_svg":
        _assert_reaction_svg(test, case, result)
    elif artifact == "connected_layout":
        _assert_connected_layout(test, result)
    elif artifact in ("png", "stereo_pair_png"):
        _assert_nonblank_image(test, result)
    elif artifact == "geometric_pair_png":
        test.assertEqual([name for name, _ in result], ["cis", "trans"])
        for _, image in result:
            _assert_nonblank_image(test, image)
    elif artifact == "docx_package":
        _assert_docx_package(test, result, len(case["input"]))
    else:
        test.fail(f"Unknown manifest artifact: {artifact}")


class PreviewManifestTests(unittest.TestCase):
    def test_manifest_cases(self):
        manifest = load_manifest()
        self.assertEqual(manifest["schema_version"], 1)
        declared_errors = {
            case["expected"]["error_type"]
            for case in manifest["cases"]
            if case["expected"]["outcome"] == "rejected"
        }
        public_errors = {
            error_type.__name__
            for error_type in structure_draw.KekuleInputError.__subclasses__()
        }
        self.assertEqual(declared_errors, public_errors)
        with tempfile.TemporaryDirectory() as temp_dir:
            for case in manifest["cases"]:
                with self.subTest(case=case["id"]):
                    expected = case["expected"]
                    if expected["outcome"] == "supported":
                        result = _call(case, temp_dir)
                        _assert_supported(self, case, result)
                        continue

                    self.assertEqual(expected["outcome"], "rejected")
                    error_type = getattr(structure_draw, expected["error_type"])
                    self.assertTrue(issubclass(error_type, structure_draw.KekuleInputError))
                    with self.assertRaises(error_type) as caught:
                        _call(case, temp_dir)
                    self.assertIn(expected["message_contains"], str(caught.exception))


if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(PreviewManifestTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if result.wasSuccessful():
        print(f"{len(load_manifest()['cases'])} preview manifest cases passed")
    sys.exit(0 if result.wasSuccessful() else 1)
