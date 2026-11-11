import json
from pathlib import Path

from dptb.entrypoints import pdso as pdso_module


def test_pdso_passes_relative_override_overlap_to_tbsystem(tmp_path, monkeypatch):
    overlap_path = tmp_path / "overlaps.h5"
    overlap_path.touch()
    structure_path = tmp_path / "structure.vasp"
    structure_path.touch()
    config_path = tmp_path / "band.json"
    config_path.write_text(json.dumps({
        "task_options": {
            "task": "band",
            "override_overlap": "./overlaps.h5",
        }
    }))

    received = {}

    class FakeTBSystem:
        def __init__(self, **kwargs):
            received.update(kwargs)

        def to_pardiso_json(self, output_dir):
            received["output_dir"] = output_dir

    monkeypatch.setattr(pdso_module, "TBSystem", FakeTBSystem)
    monkeypatch.setattr(pdso_module.subprocess, "run", lambda *args, **kwargs: None)

    pdso_module.pdso(
        INPUT=str(config_path),
        init_model="e3-model.pth",
        structure=str(structure_path),
        output_dir=str(tmp_path / "output"),
    )

    assert received["override_overlap"] == str(overlap_path.resolve())
    assert received["calculator"] == "e3-model.pth"
    assert received["data"] == str(structure_path)


def test_pdso_accepts_top_level_override_overlap(tmp_path):
    overlap_path = tmp_path / "overlaps.h5"
    config_path = tmp_path / "band.json"
    config_path.write_text(json.dumps({"override_overlap": "./overlaps.h5"}))

    resolved = pdso_module._get_override_overlap(str(config_path))

    assert resolved == str(overlap_path.resolve())
