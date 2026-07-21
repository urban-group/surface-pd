"""Tests for library diagnostics and observability boundaries."""

import ast
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
LIBRARY_MODULES_WITHOUT_STDOUT = [
    PROJECT_ROOT / "surface_pd" / "core" / "post_check.py",
    PROJECT_ROOT / "surface_pd" / "core" / "pre_check.py",
    PROJECT_ROOT / "surface_pd" / "core" / "enumeration_slab.py",
]
PACKAGE_MODULES_WITHOUT_PASS = [
    path
    for path in (PROJECT_ROOT / "surface_pd").rglob("*.py")
    if "__pycache__" not in path.parts
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


def _pass_statement_lines(path: Path) -> list[int]:
    """Return source lines that contain active pass statements."""
    tree = ast.parse(path.read_text(), filename=str(path))
    return [
        node.lineno for node in ast.walk(tree) if isinstance(node, ast.Pass)
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


def test_package_modules_do_not_use_silent_pass_statements():
    """Package control flow should make intentional no-op paths explicit."""
    pass_statements = {
        path.relative_to(PROJECT_ROOT).as_posix(): _pass_statement_lines(path)
        for path in PACKAGE_MODULES_WITHOUT_PASS
    }

    assert pass_statements == {
        path.relative_to(PROJECT_ROOT).as_posix(): []
        for path in PACKAGE_MODULES_WITHOUT_PASS
    }
