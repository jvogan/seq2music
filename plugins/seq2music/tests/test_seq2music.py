import importlib.util
import io
import json
import math
import os
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("seq2music", ROOT / "scripts" / "seq2music.py")
hj = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hj)


class Seq2MusicTests(unittest.TestCase):
    def roundtrip_meta(self, sequence, kind, step=.3):
        return {
            "profile":hj.ROUNDTRIP_VERSION,"class":"exact-normalized-sequence","mode":"encode",
            "kind":kind,"record_name":"roundtrip-test","sequence_length":len(sequence),
            "normalized_sequence_sha256":hj.sha256_bytes(sequence.encode()),"midi_channel":0,
            "ordering":"ascending-note-on-tick","scope":"normalized symbols only",
            "algorithm":hj.ALGORITHM_VERSION,"ppq":480,"step_beats":step,
            "step_ticks":max(1,round(step*480)),"gate_ticks":max(1,round(max(1,round(step*480))*.82)),
            "event_count":len(sequence),
        }

    def test_protein_pitch_is_unique_and_complete(self):
        self.assertEqual(set(hj.AA_PITCH), set(hj.AA_ORDER))
        self.assertEqual(len(set(hj.AA_PITCH.values())), 20)

    def test_slug_avoids_windows_device_names(self):
        for value in ("CON","aux.txt","LPT1","com9.mid","NUL"):
            with self.subTest(value=value):
                self.assertTrue(hj.slug(value).startswith("seq-"))

    def test_musicxml_element_budget_covers_maximum_generated_score(self):
        self.assertGreaterEqual(hj.MAX_SCORE_ELEMENTS,12*hj.MAX_MIDI_EVENTS+1024)

    @unittest.skipUnless(hasattr(os,"symlink"),"symlink support")
    def test_generated_writers_refuse_preexisting_symlink(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); victim=root/"victim"; target=root/"score.mid"; victim.write_bytes(b"SAFE")
            try: target.symlink_to(victim)
            except OSError: self.skipTest("symlink creation unavailable")
            with self.assertRaises(FileExistsError):
                hj.write_midi(target,hj.sequence_events("A","protein",.1),100)
            self.assertEqual(victim.read_bytes(),b"SAFE")

    def test_global_alignment_localizes_single_change(self):
        a, b = hj.needleman_wunsch("MKTAY", "MKAAY")
        self.assertEqual(len(a), len(b))
        self.assertEqual(sum(x != y for x, y in zip(a, b)), 1)

    def test_msa_entropy(self):
        events = hj.msa_events([("a", "AAAA"), ("b", "AACA")], "protein", 0.1)
        self.assertEqual(events[0]["entropy_bits"], 0.0)
        self.assertGreater(events[2]["entropy_bits"], 0.0)
        self.assertGreater(len(events[2]["notes"]), 1)
        self.assertNotEqual(events[2]["entropy_interval_semitones"], 12)

    def test_msa_all_gap_column_is_silent_and_not_conserved(self):
        events = hj.msa_events([("a", "A-"), ("b", "A-")], "protein", 0.1)
        self.assertEqual(events[1]["feature"], "no-informative-observation")
        self.assertEqual(events[1]["notes"], [])
        self.assertEqual(events[1]["coverage"], 0.0)

    def test_iupac_n_has_explicit_ambiguous_cue(self):
        hj.validate_sequence("ACGTN", "dna")
        event = hj.sequence_events("N", "dna", 0.1)[0]
        self.assertEqual(event["feature"], "ambiguous")
        self.assertLess(event["velocity"], 64)

    def test_structure_is_rigid_transform_invariant(self):
        residues = hj.parse_pdb(ROOT / "assets" / "examples" / "demo-structure.pdb")
        first = hj.structure_events(residues, 0.1)
        transformed = []
        for r in residues:
            # 90-degree rotation about z, followed by translation.
            transformed.append({**r, "x": -r["y"] + 11.0, "y": r["x"] - 7.0, "z": r["z"] + 3.0})
        second = hj.structure_events(transformed, 0.1)
        keep = ("notes", "velocity", "radial_rank", "contact_count")
        self.assertEqual([[e[k] for k in keep] for e in first], [[e[k] for k in keep] for e in second])

    def test_mmcif_atom_site_parser(self):
        residues = hj.parse_mmcif(ROOT / "assets" / "examples" / "demo-structure.cif")
        self.assertEqual("".join(r["symbol"] for r in residues), "MKTAYI")
        self.assertEqual({r["chain"] for r in residues}, {"A"})

    def test_mmcif_apostrophe_atom_name_does_not_break_row_tokenization(self):
        text = """data_x
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 O O5' . ALA A 1 0 0 0 1 10 1 ALA A O5' 1
ATOM 2 C CA . ALA A 1 1 2 3 1 10 1 ALA A CA 1
#
"""
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "apostrophe.cif"
            path.write_text(text)
            self.assertEqual(hj.parse_mmcif(path)[0]["symbol"], "A")

    def test_structure_neighbors_use_residue_numbers_not_observed_indices(self):
        residues = [
            {"chain":"A","resnum":1,"ins":"","symbol":"A","x":0.0,"y":0.0,"z":0.0,"bfactor":0.0},
            {"chain":"A","resnum":10,"ins":"","symbol":"A","x":1.0,"y":0.0,"z":0.0,"bfactor":0.0},
        ]
        events = hj.structure_events(residues, 0.1)
        self.assertEqual([e["contact_count"] for e in events], [1, 1])

    def test_midi_has_explicit_left_and_right_pan_controls(self):
        events = hj.diff_events("A", "C", "protein", 0.1)[0]
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "pan.mid"
            hj.write_midi(path, events, 100)
            data = path.read_bytes()
            self.assertIn(bytes([0xB0, 0x0A, 64]), data)
            self.assertIn(bytes([0xB1, 0x0A, 12]), data)
            self.assertIn(bytes([0xB2, 0x0A, 115]), data)

    def test_nucleotide_msa_n_is_ambiguity_not_fifth_base(self):
        records=[("a","A"),("c","C"),("g","G"),("t","T"),("n","N")]
        event=hj.msa_events(records,"dna",0.1)[0]
        self.assertLessEqual(event["entropy_ratio"],1.0)
        self.assertAlmostEqual(event["entropy_bits"],2.0)
        self.assertEqual(event["informative_coverage"],0.8)

    def test_mmcif_decimal_fields_are_preserved(self):
        residues=hj.parse_mmcif(ROOT / "assets" / "examples" / "demo-structure.cif")
        self.assertEqual(residues[0]["occupancy"],1.0)
        self.assertLess(max(r["bfactor"] for r in residues),1000)

    def test_structure_arbitrary_rotation_and_radial_ties(self):
        residues=[
            {"chain":"A","resnum":1,"ins":"","symbol":"A","x":-4.0,"y":0.0,"z":0.0,"bfactor":0.0},
            {"chain":"B","resnum":1,"ins":"","symbol":"A","x":4.0,"y":0.0,"z":0.0,"bfactor":0.0},
        ]
        before=hj.structure_events(residues,0.1)
        angle=math.radians(45)
        rotated=[{**r,"x":r["x"]*math.cos(angle)-r["y"]*math.sin(angle)+0.123,
                  "y":r["x"]*math.sin(angle)+r["y"]*math.cos(angle)-0.456} for r in residues]
        after=hj.structure_events(rotated,0.1)
        self.assertEqual([e["contact_count"] for e in before],[e["contact_count"] for e in after])
        self.assertEqual([e["radial_rank"] for e in before],[0.5,0.5])
        self.assertEqual([e["notes"] for e in before],[e["notes"] for e in after])

    def test_bundle_and_verify(self):
        fasta = ROOT / "assets" / "examples" / "demo-protein.fasta"
        events = hj.sequence_events(hj.read_fasta(fasta)[0][1], "protein", 0.04)
        with tempfile.TemporaryDirectory() as td:
            bundle = hj.make_bundle([fasta], events, Path(td), "sonify", "test", 100, False)
            files = {p.suffix for p in bundle.iterdir()}
            self.assertTrue({".mid", ".wav", ".csv", ".html", ".txt", ".json"} <= files)
            manifest = next(bundle.glob("*.run.json"))
            data = json.loads(manifest.read_text())
            self.assertEqual(data["algorithm"], "seq2music-mapping-v2")
            self.assertIn("verification_scope", data)
            for artifact in data["artifacts"]:
                self.assertEqual(hj.sha256_file(bundle / artifact["path"]), artifact["sha256"])
            self.assertEqual((next(bundle.glob("*.mid"))).read_bytes()[:4], b"MThd")
            self.assertEqual((next(bundle.glob("*.wav"))).read_bytes()[:4], b"RIFF")

    def test_verify_rejects_empty_and_traversal_manifests(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            empty = root / "empty.json"
            empty.write_text("{}")
            with self.assertRaises(ValueError):
                hj.cmd_verify(type("Args", (), {"manifest":str(empty)})())
            traversal = root / "traversal.json"
            traversal.write_text(json.dumps({"schema_version":"1.0","run_id":"x","artifacts":[{"path":"../outside","sha256":"0"*64}]}))
            with self.assertRaises(ValueError):
                hj.cmd_verify(type("Args", (), {"manifest":str(traversal)})())

    def test_protein_all_symbols_round_trip_from_midi_alone(self):
        sequence=hj.AA_ORDER
        with tempfile.TemporaryDirectory() as td:
            midi=Path(td)/"protein.mid"
            hj.write_midi(midi,hj.sequence_events(sequence,"protein",.18),100,self.roundtrip_meta(sequence,"protein"))
            decoded=hj.decode_midi_sequence(midi)
            self.assertEqual(decoded["sequence"],sequence)
            self.assertTrue(decoded["exact"])
            self.assertEqual(decoded["status"],"exact-embedded")
            self.assertEqual(hj.parse_midi(midi)["metadata"]["kind"],"protein")

    def test_dna_and_rna_round_trip_preserve_t_u(self):
        for kind,sequence in (("dna","ACGTN"),("rna","ACGUN")):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as td:
                midi=Path(td)/f"{kind}.mid"
                hj.write_midi(midi,hj.sequence_events(sequence,kind,.18),100,self.roundtrip_meta(sequence,kind))
                decoded=hj.decode_midi_sequence(midi)
                self.assertEqual(decoded["sequence"],sequence)
                self.assertEqual(decoded["kind"],kind)

    def test_metadata_free_midi_is_never_claimed_exact(self):
        with tempfile.TemporaryDirectory() as td:
            midi=Path(td)/"legacy.mid"
            hj.write_midi(midi,hj.sequence_events("ACGT","dna",.18),100)
            with self.assertRaises(ValueError): hj.decode_midi_sequence(midi)
            decoded=hj.decode_midi_sequence(midi,requested_kind="dna",allow_edited=True)
            self.assertEqual(decoded["sequence"],"ACGT")
            self.assertFalse(decoded["exact"])

    def test_edited_note_fails_strict_and_can_export_unverified(self):
        original="ACGT"
        with tempfile.TemporaryDirectory() as td:
            events=hj.sequence_events(original,"dna",.18); events[1]["notes"]=[hj.NUC_PITCH["G"]]
            midi=Path(td)/"edited-pitch.mid"
            hj.write_midi(midi,events,100,self.roundtrip_meta(original,"dna"))
            with self.assertRaises(ValueError): hj.decode_midi_sequence(midi)
            decoded=hj.decode_midi_sequence(midi,allow_edited=True)
            self.assertFalse(decoded["exact"])
            self.assertEqual(decoded["status"],"edited-or-unverified")

    def test_fractional_step_uses_one_shared_tick_grid(self):
        sequence="AC"
        step=.05104166875
        with tempfile.TemporaryDirectory() as td:
            midi=Path(td)/"fractional.mid"; meta=self.roundtrip_meta(sequence,"protein",step)
            hj.write_midi(midi,hj.sequence_events(sequence,"protein",step*60/100),100,meta)
            self.assertEqual([n["tick"] for n in hj.parse_midi(midi)["notes"]],[0,meta["step_ticks"]])
            self.assertEqual(hj.decode_midi_sequence(midi)["sequence"],sequence)

    def test_manifest_matching_upgrades_status_and_mismatch_fails(self):
        sequence="MKTAY"
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); midi=root/"sample.mid"
            meta=self.roundtrip_meta(sequence,"protein")
            hj.write_midi(midi,hj.sequence_events(sequence,"protein",.18),100,meta)
            manifest=root/"sample.run.json"
            manifest.write_text(json.dumps({"roundtrip":meta,"artifacts":[{"path":midi.name,"sha256":hj.sha256_file(midi)}]}))
            self.assertEqual(hj.decode_midi_sequence(midi,manifest)["status"],"exact-manifest-matched")
            data=json.loads(manifest.read_text()); data["artifacts"][0]["sha256"]="0"*64; manifest.write_text(json.dumps(data))
            with self.assertRaises(ValueError): hj.decode_midi_sequence(midi,manifest)
            self.assertFalse(hj.decode_midi_sequence(midi,manifest,allow_edited=True)["exact"])

    def test_encode_decode_cli_writes_fixed_safe_outputs(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); source=root/"sample.fasta"; source.write_text(">sample\nMKTAYIAK\n")
            rendered=root/"rendered"
            with redirect_stdout(io.StringIO()):
                self.assertEqual(hj.main(["encode","--input",str(source),"--kind","protein","--out",str(rendered)]),0)
            bundle=next(rendered.iterdir()); midi=next(bundle.glob("*.mid")); manifest=next(bundle.glob("*.run.json")); decoded=root/"decoded"
            with redirect_stdout(io.StringIO()):
                self.assertEqual(hj.main(["decode","--midi",str(midi),"--manifest",str(manifest),"--out",str(decoded)]),0)
            self.assertIn("MKTAYIAK",(decoded/"sequence.fasta").read_text())
            report=json.loads((decoded/"decode.report.json").read_text())
            self.assertTrue(report["exact"])
            self.assertEqual(report["status"],"exact-manifest-matched")
            with redirect_stdout(io.StringIO()):
                self.assertEqual(hj.main(["decode","--midi",str(midi),"--out",str(decoded)]),2)

    def test_msa_and_structure_midi_are_not_decodable_as_exact_sequences(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for name,events in (
                ("msa",hj.msa_events([("a","AA"),("b","AC")],"protein",.18)),
                ("structure",hj.structure_events(hj.parse_pdb(ROOT/"assets/examples/demo-structure.pdb"),.18)),
            ):
                midi=root/f"{name}.mid"; hj.write_midi(midi,events,100)
                with self.assertRaises(ValueError): hj.decode_midi_sequence(midi)

    def test_headerless_long_filename_is_bounded_and_roundtrips(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); source=root/("x"*210+".fa"); source.write_text("MKTAY\n")
            records=hj.read_fasta(source)
            self.assertEqual(len(records[0][0]),200)
            with redirect_stdout(io.StringIO()):
                self.assertEqual(hj.main(["encode","--input",str(source),"--kind","protein","--out",str(root/"out")]),0)
            bundle=next((root/"out").iterdir())
            self.assertLessEqual(len(bundle.name.split("-encode-")[0]),96)
            decoded=hj.decode_midi_sequence(next(bundle.glob("*.mid")))
            self.assertEqual(decoded["sequence"],"MKTAY")

    def test_extended_iupac_is_msa_only_and_not_counted_as_bases(self):
        records=[("a","A"),("c","C"),("r","R"),("y","Y"),("gap","-")]
        event=hj.msa_events(records,"dna",.1)[0]
        self.assertEqual(event["entropy_bits"],1.0)
        self.assertEqual(event["entropy_ratio"],.5)
        self.assertEqual(event["coverage"],.8)
        self.assertEqual(event["informative_coverage"],.4)
        self.assertEqual(event["ambiguous_coverage"],.4)
        silent=hj.msa_events([("r","R"),("y","Y"),("n","N")],"dna",.1)[0]
        self.assertEqual(silent["notes"],[])
        self.assertEqual(silent["symbol"],"?")
        self.assertEqual(silent["ambiguous_symbols"],"NRY")
        voiced=hj.msa_events([("a","A"),("r","R"),("y","Y")],"dna",.1)[0]
        self.assertEqual(voiced["ambiguous_symbols"],"RY")
        self.assertEqual(voiced["ambiguous_count"],2)
        hj.validate_sequence("ACGTNRYSWKMBDHV-","dna",allow_gaps=True,allow_iupac=True)
        hj.validate_sequence("ACGUNRYSWKMBDHV-","rna",allow_gaps=True,allow_iupac=True)
        with self.assertRaises(ValueError): hj.validate_sequence("ACGR","dna")
        with self.assertRaises(ValueError): hj.validate_sequence("ACGU","dna",allow_iupac=True)
        with self.assertRaises(ValueError): hj.validate_sequence("ACGT","rna",allow_iupac=True)
        self.assertEqual(hj.infer_kind("MKTAR","auto"),"protein")
        with self.assertRaises(ValueError): hj.msa_events([("a","A?"),("b","AA")],"dna",.1)

    def test_diff_preflights_exact_event_and_matrix_limits(self):
        with self.assertRaisesRegex(ValueError,"1,002,001"):
            hj.needleman_wunsch("A"*1000,"A"*1000)
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); a=root/"a.fa"; b=root/"b.fa"
            a.write_text(">a\nAAAA\n"); b.write_text(">b\nAAA\n")
            stderr=io.StringIO()
            with redirect_stdout(io.StringIO()), redirect_stderr(stderr):
                self.assertEqual(hj.main(["diff","--reference",str(a),"--variant",str(b),"--kind","protein","--max-events","6","--out",str(root/"out")]),2)
            self.assertIn("7 events",stderr.getvalue())

    def test_structure_model_selection_for_pdb_and_mmcif(self):
        def pdb_atom(model_x):
            return f"ATOM      1  CA  ALA A   1    {model_x:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{10:6.2f}           C  "
        pdb_text="MODEL        1\n"+pdb_atom(1)+"\nENDMDL\nMODEL        2\n"+pdb_atom(9)+"\nENDMDL\n"
        cif_text="""data_x
loop_
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
CA ALA A 1 1 0 0 1 10 1
CA ALA A 1 9 0 0 1 10 2
#
"""
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); pdb=root/"two.pdb"; cif=root/"two.cif"
            pdb.write_text(pdb_text); cif.write_text(cif_text)
            self.assertEqual(hj.parse_pdb(pdb)[0]["model"],1)
            self.assertEqual(hj.parse_pdb(pdb,model_selection=2)[0]["x"],9.0)
            self.assertEqual(hj.parse_mmcif(cif)[0]["model"],1)
            self.assertEqual(hj.parse_mmcif(cif,model_selection=2)[0]["x"],9.0)
            implicit=root/"implicit.pdb"; implicit.write_text(pdb_atom(1)+"\n")
            self.assertEqual(hj.parse_pdb(implicit,model_selection=1)[0]["model"],1)
            with self.assertRaisesRegex(ValueError,"requested PDB model 2 is not present"):
                hj.parse_pdb(implicit,model_selection=2)
            mixed=root/"mixed.pdb"; mixed.write_text(pdb_atom(1)+"\nMODEL        2\n"+pdb_atom(9)+"\nENDMDL\n")
            with self.assertRaisesRegex(ValueError,"outside a MODEL"):
                hj.parse_pdb(mixed)

    def test_mmcif_model_numbers_fail_explicitly_and_do_not_fallback(self):
        template="""data_x
loop_
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
{atom} ALA A 1 0 0 0 {model}
#
"""
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"model.cif"
            for bad in ("1.9","NaN","inf","0"):
                path.write_text(template.format(atom="CA",model=bad))
                with self.subTest(model=bad), self.assertRaisesRegex(ValueError,"positive integer"):
                    hj.parse_mmcif(path)
            for missing in (".","?"):
                path.write_text(template.format(atom="CA",model=missing))
                with self.subTest(model=missing), self.assertRaisesRegex(ValueError,"is missing"):
                    hj.parse_mmcif(path)
            path.write_text(template.format(atom="CB",model="2"))
            with self.assertRaisesRegex(ValueError,"no supported protein C-alpha atoms.*requested model 2"):
                hj.parse_mmcif(path,model_selection=2)

    def test_implicit_mmcif_model_one_and_pdb_present_model_without_ca(self):
        cif="""data_x
loop_
_atom_site.auth_atom_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
CA ALA A 1 0 0 0
#
"""
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); path=root/"implicit.cif"; path.write_text(cif)
            self.assertEqual(hj.parse_mmcif(path,model_selection=1)[0]["model"],1)
            with self.assertRaisesRegex(ValueError,"only implicit model 1"):
                hj.parse_mmcif(path,model_selection=2)
            pdb=root/"no-ca.pdb"; pdb.write_text("MODEL        2\nATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 10.00           C  \nENDMDL\n")
            with self.assertRaisesRegex(ValueError,"no supported protein C-alpha atoms.*requested model 2"):
                hj.parse_pdb(pdb,model_selection=2)

    def test_structure_cli_records_selected_model(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); source=ROOT/"assets/examples/demo-structure.pdb"
            with redirect_stdout(io.StringIO()):
                self.assertEqual(hj.main(["structure","--input",str(source),"--model","1","--out",str(root)]),0)
            bundle=next(root.iterdir()); manifest=json.loads(next(bundle.glob("*.run.json")).read_text())
            self.assertEqual(manifest["details"]["model"],1)
            self.assertEqual(manifest["parameters"]["selected_model"],1)
            self.assertIn("model",next(bundle.glob("*.events.csv")).read_text().splitlines()[0])

    def test_manifests_use_strict_json(self):
        sequence="ACGT"
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); midi=root/"sample.mid"; meta=self.roundtrip_meta(sequence,"dna")
            hj.write_midi(midi,hj.sequence_events(sequence,"dna",.1),100,meta)
            duplicate=root/"duplicate.json"
            duplicate.write_text('{"roundtrip":{},"roundtrip":{},"artifacts":[]}')
            with self.assertRaisesRegex(ValueError,"duplicate JSON key"):
                hj.decode_midi_sequence(midi,duplicate)
            nan=root/"nan.json"
            nan.write_text('{"schema_version":"1.0","run_id":"x","dead":NaN,"artifacts":[]}')
            with self.assertRaisesRegex(ValueError,"invalid JSON constant"):
                hj.cmd_verify(type("Args",(),{"manifest":str(nan)})())
            null_artifacts=root/"null-artifacts.json"
            null_artifacts.write_text(json.dumps({"roundtrip":meta,"artifacts":None}))
            with self.assertRaisesRegex(ValueError,"artifacts must be a list"):
                hj.decode_midi_sequence(midi,null_artifacts)
            bad_digest=root/"bad-digest.json"
            bad_digest.write_text(json.dumps({"roundtrip":meta,"artifacts":[{"path":midi.name,"sha256":"not-a-digest"}]}))
            with self.assertRaisesRegex(ValueError,"invalid artifact"):
                hj.decode_midi_sequence(midi,bad_digest,allow_edited=True)

    def test_roundtrip_magic_is_only_accepted_in_meta_type_7f(self):
        sequence="AC"
        with tempfile.TemporaryDirectory() as td:
            midi=Path(td)/"wrong-meta.mid"
            hj.write_midi(midi,hj.sequence_events(sequence,"protein",.1),100,self.roundtrip_meta(sequence,"protein"))
            data=midi.read_bytes(); marker=b"\xff\x7f"+hj.varlen(len(hj.ROUNDTRIP_MAGIC+json.dumps(self.roundtrip_meta(sequence,"protein"),sort_keys=True,separators=(",",":"),ensure_ascii=True).encode()))
            self.assertIn(marker,data)
            midi.write_bytes(data.replace(marker,b"\xff\x01"+marker[2:],1))
            with self.assertRaisesRegex(ValueError,"meta type 0x7F"):
                hj.parse_midi(midi)

    def test_decode_hashes_the_midi_snapshot_once(self):
        sequence="AC"
        with tempfile.TemporaryDirectory() as td:
            midi=Path(td)/"once.mid"
            hj.write_midi(midi,hj.sequence_events(sequence,"protein",.1),100,self.roundtrip_meta(sequence,"protein"))
            original=hj.sha256_file
            def reject_midi_reopen(path,*args,**kwargs):
                if Path(path)==midi: raise AssertionError("decoder reopened MIDI for hashing")
                return original(path,*args,**kwargs)
            hj.sha256_file=reject_midi_reopen
            try: self.assertEqual(hj.decode_midi_sequence(midi)["sequence"],sequence)
            finally: hj.sha256_file=original

    @unittest.skipUnless(hj.supports_anchored_output_transaction(),"secure forced replacement requires directory descriptors")
    def test_force_replaces_only_known_generated_directory(self):
        fasta=ROOT/"assets/examples/demo-protein.fasta"
        events=hj.sequence_events(hj.read_fasta(fasta)[0][1],"protein",.04)
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); bundle=hj.make_bundle([fasta],events,root,"sonify","test",100,False)
            (bundle/"user-notes.txt").write_text("keep")
            with self.assertRaisesRegex(ValueError,"unknown or unsafe"):
                hj.make_bundle([fasta],events,root,"sonify","test",100,True)
            self.assertEqual((bundle/"user-notes.txt").read_text(),"keep")

    @unittest.skipUnless(hasattr(os,"geteuid"),"POSIX ownership check")
    def test_output_parent_rejects_shared_writable_ancestry(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); shared=root/"shared"; child=shared/"owned"
            shared.mkdir(); child.mkdir(); shared.chmod(0o777)
            try:
                with self.assertRaisesRegex(ValueError,"writable by another account"):
                    hj.ensure_real_directory(child)
            finally:
                shared.chmod(0o700)

    def test_musicxml_uses_one_bounded_streaming_parse(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); score=root/"score.musicxml"
            score.write_text("<score-partwise><part id='P1'><measure><note><pitch><step>C</step><octave>4</octave></pitch></note></measure></part></score-partwise>")
            original=hj.ET.fromstring
            hj.ET.fromstring=lambda *_args,**_kwargs: (_ for _ in ()).throw(AssertionError("full-tree reparse"))
            try: self.assertEqual(hj.parse_musicxml(score)["onsets"],[[60]])
            finally: hj.ET.fromstring=original
            original_limit=hj.MAX_SCORE_ELEMENTS; hj.MAX_SCORE_ELEMENTS=8
            score.write_text("<score-partwise><part id='P1'><measure><credit/><credit/><credit/><credit/><note><pitch><step>C</step><octave>4</octave></pitch></note></measure></part></score-partwise>")
            try:
                with self.assertRaisesRegex(ValueError,"too many XML elements"):
                    hj.parse_musicxml(score)
            finally: hj.MAX_SCORE_ELEMENTS=original_limit

    def test_portable_publication_is_new_output_only(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); staged=root/"staged"; final=root/"final"; staged.mkdir(); (staged/"known.txt").write_text("new")
            original=hj.supports_anchored_output_transaction; hj.supports_anchored_output_transaction=lambda:False
            try:
                hj.commit_generated_directory(staged,final,{"known.txt"},False)
                self.assertEqual((final/"known.txt").read_text(),"new")
                staged=root/"staged-again"; staged.mkdir(); (staged/"known.txt").write_text("replacement")
                with self.assertRaisesRegex(ValueError,"requires secure directory-descriptor support"):
                    hj.commit_generated_directory(staged,final,{"known.txt"},True)
                self.assertEqual((final/"known.txt").read_text(),"new")
                self.assertEqual((staged/"known.txt").read_text(),"replacement")
            finally: hj.supports_anchored_output_transaction=original

    def test_portable_publication_rejects_unexpected_entries(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); staged=root/"staged"; final=root/"final"; staged.mkdir()
            (staged/"known.txt").write_text("new"); (staged/"unexpected.txt").write_text("no")
            original=hj.supports_anchored_output_transaction; hj.supports_anchored_output_transaction=lambda:False
            try:
                with self.assertRaisesRegex(ValueError,"unknown or unsafe"):
                    hj.commit_generated_directory(staged,final,{"known.txt"},False)
            finally: hj.supports_anchored_output_transaction=original

    @unittest.skipUnless(hj.supports_anchored_output_transaction(),"directory-descriptor transaction test")
    def test_force_rechecks_directory_after_identity_anchored_move(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); staged=root/"staged"; final=root/"final"
            staged.mkdir(); final.mkdir(); (staged/"known.txt").write_text("new"); (final/"known.txt").write_text("old")
            original=hj.os.rename; injected=False
            def inject_before_move(src,dst,*args,**kwargs):
                nonlocal injected
                if src==final.name and not injected:
                    (final/"DO-NOT-DELETE").write_text("user")
                    injected=True
                return original(src,dst,*args,**kwargs)
            hj.os.rename=inject_before_move
            try:
                with self.assertRaisesRegex(ValueError,"unknown or unsafe"):
                    hj.commit_generated_directory(staged,final,{"known.txt"},True)
            finally: hj.os.rename=original
            self.assertEqual((final/"DO-NOT-DELETE").read_text(),"user")
            self.assertEqual((final/"known.txt").read_text(),"old")
            self.assertTrue(staged.exists())

    @unittest.skipUnless(hj.supports_anchored_output_transaction(),"directory-descriptor transaction test")
    def test_swapped_staging_directory_is_quarantined_and_old_output_restored(self):
        for has_existing in (False,True):
            with self.subTest(has_existing=has_existing), tempfile.TemporaryDirectory() as td:
                root=Path(td); staged=root/"staged"; final=root/"final"; staged.mkdir(); (staged/"known.txt").write_text("legitimate")
                if has_existing:
                    final.mkdir(); (final/"known.txt").write_text("old")
                original=hj.os.rename; swapped=False
                def swap_stage(src,dst,*args,**kwargs):
                    nonlocal swapped
                    if src==staged.name and dst==final.name and not swapped:
                        swapped=True
                        original(staged.name,"legitimate-staged",*args,**kwargs)
                        staged.mkdir(); (staged/"known.txt").write_text("attacker")
                    return original(src,dst,*args,**kwargs)
                hj.os.rename=swap_stage
                try:
                    with self.assertRaisesRegex(ValueError,"quarantined"):
                        hj.commit_generated_directory(staged,final,{"known.txt"},has_existing)
                finally: hj.os.rename=original
                if has_existing:
                    self.assertEqual((final/"known.txt").read_text(),"old")
                else:
                    self.assertFalse(final.exists())
                quarantines=list(root.glob(".final.quarantine-*"))
                self.assertEqual(len(quarantines),1)
                self.assertEqual((quarantines[0]/"known.txt").read_text(),"attacker")

    def test_long_explicit_decode_output_uses_bounded_temp_prefix(self):
        sequence="AC"
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); midi=root/"sample.mid"
            hj.write_midi(midi,hj.sequence_events(sequence,"protein",.1),100,self.roundtrip_meta(sequence,"protein"))
            out=root/("d"*180)
            with redirect_stdout(io.StringIO()):
                self.assertEqual(hj.main(["decode","--midi",str(midi),"--out",str(out)]),0)
            self.assertTrue((out/"sequence.fasta").is_file())

    def test_generated_wav_and_musicxml_roundtrip_all_kinds(self):
        for kind,sequence in (("protein","ACDEFG"),("dna","ACGTN"),("rna","ACGUN")):
            with self.subTest(kind=kind), tempfile.TemporaryDirectory() as td:
                root=Path(td); source=root/f"{kind}.fasta"; source.write_text(f">{kind}\n{sequence}\n")
                with redirect_stdout(io.StringIO()):
                    self.assertEqual(hj.main(["encode","--input",str(source),"--kind",kind,"--out",str(root/"rendered")]),0)
                bundle=next((root/"rendered").iterdir()); wav=next(bundle.glob("*.wav")); score=next(bundle.glob("*.musicxml")); svg=next(bundle.glob("*.score.svg"))
                audio=hj.decode_wave_sequence(wav); notation=hj.decode_musicxml_sequence(score)
                self.assertEqual((audio["sequence"],audio["status"],audio["exact"]),(sequence,"exact-wav-embedded",True))
                self.assertEqual((notation["sequence"],notation["status"],notation["exact"]),(sequence,"exact-score-embedded",True))
                self.assertEqual(audio["kind_source"],"embedded"); self.assertEqual(notation["kind_source"],"embedded")
                self.assertTrue(audio["confidence"]["not_applicable"]); self.assertTrue(notation["confidence"]["not_applicable"])
                svg_text=svg.read_text()
                self.assertIn("<title>",svg_text); self.assertIn("<desc>",svg_text); self.assertIn('class="staff"',svg_text)
                manifest=json.loads(next(bundle.glob("*.run.json")).read_text())
                artifact_names={a["path"] for a in manifest["artifacts"]}
                self.assertIn(wav.name,artifact_names); self.assertIn(score.name,artifact_names); self.assertIn(svg.name,artifact_names)
                self.assertNotIn("normalized_sequence",manifest["roundtrip"])
                self.assertNotIn(b'"normalized_sequence"',next(bundle.glob("*.mid")).read_bytes())
                for measure in [e for e in hj.ET.parse(score).getroot().iter() if hj.xml_local(e)=="measure"]:
                    durations=[int(next(c for c in note if hj.xml_local(c)=="duration").text) for note in measure
                               if hj.xml_local(note)=="note" and not any(hj.xml_local(c)=="chord" for c in note)]
                    self.assertLessEqual(sum(durations),4)

    def test_plain_wav_quantizes_to_default_protein_and_visual_output(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); wav=root/"plain.wav"; expected="ACDE"
            hj.write_wav(wav,hj.sequence_events(expected,"protein",.24))
            decoded=hj.decode_wave_sequence(wav,window_seconds=.08,silence_threshold=.005)
            self.assertEqual(decoded["sequence"],expected)
            self.assertEqual(decoded["kind"],"protein")
            self.assertEqual(decoded["kind_source"],"default-protein")
            self.assertEqual(decoded["status"],"lossy-audio-quantized")
            self.assertEqual(decoded["source"]["analysis"]["window_seconds"],.08)
            with redirect_stdout(io.StringIO()):
                self.assertEqual(hj.main(["audio-decode","--audio",str(wav),"--window-seconds","0.08","--silence-threshold","0.005","--out",str(root/"decoded")]),0)
            self.assertTrue((root/"decoded/score.svg").is_file())
            self.assertFalse(json.loads((root/"decoded/decode.report.json").read_text())["exact"])

    def test_wav_pcm_tamper_requires_allow_edited(self):
        sequence="ACDE"
        with tempfile.TemporaryDirectory() as td:
            wav=Path(td)/"encoded.wav"; meta=self.roundtrip_meta(sequence,"protein"); meta["normalized_sequence"]=sequence
            hj.write_wav(wav,hj.sequence_events(sequence,"protein",.24),roundtrip=meta)
            data=bytearray(wav.read_bytes()); pos=data.find(b"data")
            self.assertGreater(pos,0); payload=pos+8; data[payload+20]^=1; wav.write_bytes(data)
            with self.assertRaisesRegex(ValueError,"PCM does not match"):
                hj.decode_wave_sequence(wav)
            decoded=hj.decode_wave_sequence(wav,allow_edited=True,window_seconds=.08,silence_threshold=.005)
            self.assertFalse(decoded["exact"]); self.assertEqual(decoded["status"],"lossy-audio-quantized")

    def test_pcm_widths_and_antiphase_stereo_are_supported(self):
        sample_rate=8000; count=1600
        for width in (1,2,3,4):
            with self.subTest(width=width), tempfile.TemporaryDirectory() as td:
                path=Path(td)/"tone.wav"; frames=bytearray()
                for index in range(count):
                    value=math.sin(2*math.pi*440*index/sample_rate)*.6
                    if width==1: frames.append(max(0,min(255,round(128+value*127))))
                    else:
                        integer=round(value*((1<<(width*8-1))-1)); frames.extend(int(integer).to_bytes(width,"little",signed=True))
                with hj.wave.open(str(path),"wb") as out:
                    out.setnchannels(1); out.setsampwidth(width); out.setframerate(sample_rate); out.writeframes(frames)
                self.assertEqual(hj.decode_wave_sequence(path,window_seconds=.08,silence_threshold=.005)["sequence"],"P")
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"antiphase.wav"; frames=bytearray()
            for index in range(count):
                value=round(math.sin(2*math.pi*440*index/sample_rate)*20000)
                frames.extend(hj.struct.pack("<hh",value,-value))
            with hj.wave.open(str(path),"wb") as out:
                out.setnchannels(2); out.setsampwidth(2); out.setframerate(sample_rate); out.writeframes(frames)
            self.assertEqual(hj.decode_wave_sequence(path,window_seconds=.08,silence_threshold=.005)["sequence"],"P")

    def test_zero_threshold_does_not_turn_digital_silence_into_a_symbol(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"silence.wav"
            with hj.wave.open(str(path),"wb") as out:
                out.setnchannels(1); out.setsampwidth(2); out.setframerate(8000); out.writeframes(b"\x00\x00"*800)
            with self.assertRaisesRegex(ValueError,"no voiced windows"):
                hj.decode_wave_sequence(path,window_seconds=.08,silence_threshold=0)

    def test_wav_duration_is_not_artificially_capped_at_sixty_seconds(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"long.wav"; sequence="A"; meta=self.roundtrip_meta(sequence,"protein"); meta["normalized_sequence"]=sequence
            event={**hj.sequence_events(sequence,"protein",61.0)[0],"duration":61.0}
            hj.write_wav(path,[event],sample_rate=8000,roundtrip=meta,max_bytes=8*hj.MIB)
            parsed=hj.parse_wave(path,max_bytes=8*hj.MIB)
            try: self.assertGreater(parsed["duration_seconds"],60)
            finally: hj.close_wave_snapshot(parsed)
            decoded=hj.decode_wave_sequence(path,max_input_bytes=8*hj.MIB)
            self.assertTrue(decoded["exact"]); self.assertEqual(decoded["sequence"],sequence)

    def test_wav_byte_budget_is_explicit_and_adjustable(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"budget.wav"
            event={**hj.sequence_events("A","protein",61.0)[0],"duration":61.0}
            with self.assertRaisesRegex(ValueError,"--max-wav-mib"):
                hj.write_wav(path,[event],sample_rate=8000,max_bytes=hj.MIB)

    def test_score_edit_fails_strict_and_allows_lossy_transcription(self):
        sequence="AC"
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); score=root/"score.musicxml"; meta=self.roundtrip_meta(sequence,"protein"); meta["normalized_sequence"]=sequence
            hj.write_musicxml(score,"test",hj.sequence_events(sequence,"protein",.2),meta)
            score.write_text(score.read_text().replace("<step>C</step>","<step>B</step>",1))
            with self.assertRaisesRegex(ValueError,"do not match"):
                hj.decode_musicxml_sequence(score)
            decoded=hj.decode_musicxml_sequence(score,allow_edited=True)
            self.assertFalse(decoded["exact"]); self.assertEqual(decoded["status"],"lossy-score-transcription")

    def test_score_visible_labels_and_timing_are_part_of_strict_semantics(self):
        sequence="T"
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); score=root/"score.musicxml"; meta=self.roundtrip_meta(sequence,"dna"); meta["normalized_sequence"]=sequence
            hj.write_musicxml(score,"test",hj.sequence_events(sequence,"dna",.2),meta)
            score.write_text(score.read_text().replace("<text>T</text>","<text>U</text>"))
            with self.assertRaisesRegex(ValueError,"do not match"):
                hj.decode_musicxml_sequence(score)
            self.assertFalse(hj.decode_musicxml_sequence(score,allow_edited=True)["exact"])
            score.write_text(score.read_text().replace("<text>U</text>","<text>T</text>").replace("<note>","<backup><duration>1</duration></backup><note>",1))
            with self.assertRaisesRegex(ValueError,"backup/forward"):
                hj.decode_musicxml_sequence(score,allow_edited=True)

    def test_exact_wave_and_score_require_embedded_kind(self):
        sequence="A"
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); meta=self.roundtrip_meta(sequence,"protein"); meta["normalized_sequence"]=sequence; meta.pop("kind")
            wav=root/"bad.wav"; score=root/"bad.musicxml"; events=hj.sequence_events(sequence,"protein",.2)
            hj.write_wav(wav,events,roundtrip=meta); hj.write_musicxml(score,"test",events,meta)
            with self.assertRaisesRegex(ValueError,"molecule kind"):
                hj.decode_wave_sequence(wav)
            with self.assertRaisesRegex(ValueError,"molecule kind"):
                hj.decode_musicxml_sequence(score)

    def test_import_receipts_enforce_the_event_safety_limit(self):
        original_limit=hj.MAX_EVENTS_HARD; hj.MAX_EVENTS_HARD=10
        try:
            sequence="A"*(hj.MAX_EVENTS_HARD+1)
            with tempfile.TemporaryDirectory() as td:
                root=Path(td); meta=self.roundtrip_meta(sequence,"protein"); meta["normalized_sequence"]=sequence
                wav=root/"too-many.wav"
                hj.write_wav(wav,hj.sequence_events("A","protein",.01),roundtrip=meta)
                with self.assertRaisesRegex(ValueError,"missing required fields|event safety limit"):
                    hj.decode_wave_sequence(wav)
                notes="".join("<note><pitch><step>C</step><octave>4</octave></pitch></note>" for _ in range(hj.MAX_EVENTS_HARD+1))
                score=root/"too-many.musicxml"
                score.write_text(f'<score-partwise><part-list><score-part id="P1"><part-name>x</part-name></score-part></part-list><part id="P1"><measure>{notes}</measure></part></score-partwise>')
                with self.assertRaisesRegex(ValueError,"event safety limit"):
                    hj.decode_musicxml_sequence(score)
        finally: hj.MAX_EVENTS_HARD=original_limit

    def test_metadata_free_musicxml_defaults_protein_and_chords_need_policy(self):
        score_text="""<?xml version="1.0"?><score-partwise version="4.0"><part-list><score-part id="P1"><part-name>x</part-name></score-part></part-list><part id="P1"><measure number="1"><note><pitch><step>C</step><octave>4</octave></pitch></note><note><chord/><pitch><step>E</step><octave>4</octave></pitch></note></measure></part></score-partwise>"""
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"chord.musicxml"; path.write_text(score_text)
            with self.assertRaisesRegex(ValueError,"contains chords"):
                hj.decode_musicxml_sequence(path)
            decoded=hj.decode_musicxml_sequence(path,polyphony="lowest")
            self.assertEqual(decoded["kind"],"protein"); self.assertEqual(decoded["sequence"],"I")
            self.assertEqual(decoded["status"],"lossy-score-transcription")

    def test_musicxml_ties_are_rejected_instead_of_double_counted(self):
        score_text='''<score-partwise><part-list><score-part id="P1"><part-name>x</part-name></score-part></part-list><part id="P1"><measure><note><pitch><step>C</step><octave>4</octave></pitch><tie type="start"/></note><note><pitch><step>C</step><octave>4</octave></pitch><tie type="stop"/></note></measure></part></score-partwise>'''
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"tie.musicxml"; path.write_text(score_text)
            with self.assertRaisesRegex(ValueError,"note semantics"):
                hj.decode_musicxml_sequence(path)

    def test_musicxml_rejects_dtd_and_duplicate_codec_metadata(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); dtd=root/"dtd.musicxml"; dtd.write_text('<!DOCTYPE x [<!ENTITY y "z">]><score-partwise/>')
            with self.assertRaisesRegex(ValueError,"DTD/entity"):
                hj.parse_musicxml(dtd)
            utf16=root/"utf16.musicxml"; utf16.write_bytes('<!DOCTYPE x [<!ENTITY y "z">]><score-partwise/>'.encode("utf-16"))
            with self.assertRaisesRegex(ValueError,"UTF-16"):
                hj.parse_musicxml(utf16)
            sequence="AC"; score=root/"duplicate.musicxml"; meta=self.roundtrip_meta(sequence,"protein"); meta["normalized_sequence"]=sequence
            hj.write_musicxml(score,"test",hj.sequence_events(sequence,"protein",.2),meta)
            tree=hj.ET.parse(score); misc=next(e for e in tree.getroot().iter() if hj.xml_local(e)=="miscellaneous")
            original=next(iter(misc)); misc.append(hj.ET.fromstring(hj.ET.tostring(original)))
            tree.write(score,encoding="utf-8",xml_declaration=True)
            with self.assertRaisesRegex(ValueError,"at most one"):
                hj.parse_musicxml(score)

    def test_musicxml_timewise_counts_one_logical_part_across_measures(self):
        score='''<score-timewise><part-list><score-part id="P1"><part-name>x</part-name></score-part></part-list><measure><part id="P1"><note><pitch><step>C</step><octave>4</octave></pitch></note></part></measure><measure><part id="P1"><note><pitch><step>D</step><octave>4</octave></pitch></note></part></measure></score-timewise>'''
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"timewise.musicxml"; path.write_text(score)
            decoded=hj.decode_musicxml_sequence(path)
            self.assertEqual(decoded["sequence"],"IK")
            self.assertEqual(decoded["source"]["parts"],1)

    def test_codec_metadata_rejects_non_string_name_and_bool_count(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); sequence="A"; meta=self.roundtrip_meta(sequence,"protein"); meta["normalized_sequence"]=sequence; meta["record_name"]=17
            wav=root/"bad-name.wav"; score=root/"bad-count.musicxml"; events=hj.sequence_events(sequence,"protein",.2)
            hj.write_wav(wav,events,roundtrip=meta)
            with self.assertRaisesRegex(ValueError,"record_name"):
                hj.decode_wave_sequence(wav)
            meta["record_name"]="ok"; hj.write_musicxml(score,"test",events,meta)
            tree=hj.ET.parse(score); field=next(e for e in tree.getroot().iter() if hj.xml_local(e)=="miscellaneous-field")
            payload=json.loads(hj.base64.b64decode(field.text)); payload["sequence_length"]=True
            field.text=hj.base64.b64encode(json.dumps(payload,separators=(",",":")).encode()).decode()
            tree.write(score,encoding="utf-8",xml_declaration=True)
            with self.assertRaisesRegex(ValueError,"missing required fields"):
                hj.decode_musicxml_sequence(score)


if __name__ == "__main__":
    unittest.main()
