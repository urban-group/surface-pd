"""Tests for library diagnostics and observability boundaries."""

import ast
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LIBRARY_MODULES_WITHOUT_STDOUT = [
    PROJECT_ROOT / "surface_pd" / "core" / "post_check.py",
    PROJECT_ROOT / "surface_pd" / "core" / "pre_check.py",
    PROJECT_ROOT / "surface_pd" / "core" / "slab.py",
    PROJECT_ROOT / "surface_pd" / "plot" / "pd_data.py",
]


def _print_call_lines(path: Path) -> list[int]:
    """Return source lines that call the builtin print function."""
    tree = ast.parse(path.read_text(), filename=str(path))
    return [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "print"
    ]


def test_library_modules_do_not_print_to_stdout():
    """Library diagnostics should use logging or exceptions, not stdout."""
    print_calls = {
        path.relative_to(PROJECT_ROOT).as_posix(): _print_call_lines(path)
        for path in LIBRARY_MODULES_WITHOUT_STDOUT
    }

    assert print_calls == {
        path.relative_to(PROJECT_ROOT).as_posix(): []
        for path in LIBRARY_MODULES_WITHOUT_STDOUT
    }
