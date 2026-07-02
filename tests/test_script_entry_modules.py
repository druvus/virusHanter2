"""Guard: Snakemake ``script:`` entry modules must not use a future import.

Snakemake prepends its own boilerplate when materialising a ``script:``
rule, so a ``from __future__ import annotations`` line no longer sits at
the top of the file and the interpreter raises

    SyntaxError: from __future__ imports must occur at the beginning of
    the file

This only surfaces at pipeline runtime (direct ``import`` of the module
is unaffected), so a unit-level guard is the cheapest way to catch a
regression. The affected modules are the ones a rule invokes via
``script:``.
"""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"

# Modules invoked by a Snakemake `script:` directive.
SCRIPT_ENTRY_MODULES = [
    "collect_software_versions.py",
    "write_provenance.py",
    "aggregate_run_information.py",
]


@pytest.mark.parametrize("module", SCRIPT_ENTRY_MODULES)
def test_script_entry_module_has_no_future_import(module):
    # Parse the AST so a docstring/comment that merely *mentions* the
    # future import (as these modules do, to explain the constraint) is
    # not mistaken for an actual statement.
    tree = ast.parse((SCRIPTS / module).read_text())
    future = [
        node
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module == "__future__"
    ]
    assert not future, (
        f"{module} is run via Snakemake script:; a __future__ import breaks "
        "at runtime because Snakemake prepends its own preamble."
    )
