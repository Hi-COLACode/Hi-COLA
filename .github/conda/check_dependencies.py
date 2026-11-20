"""Check the conda recipe's run requirements match pyproject.toml.

Every Python dependency in pyproject.toml must appear in meta.yaml's run
requirements with the same version specifier. Other run requirements must be
python or a native library also listed under host.
"""

import logging
from pathlib import Path
import re
import sys
import tomllib

logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parents[2]
# PyPI names whose conda-forge package has a different name.
CONDA_NAMES = {"matplotlib": "matplotlib-base"}
REQUIREMENT = re.compile(r"\s*([A-Za-z0-9_.\-]+)\s*(.*?)\s*$")


def split(requirement):
    """Return (name, specifier) with whitespace removed from the specifier."""
    name, specifier = REQUIREMENT.match(requirement).groups()
    return name.lower(), specifier.replace(" ", "")


def recipe_section(text, section):
    """Return the requirements listed under requirements/<section> in meta.yaml."""
    entries, inside = [], False
    for line in text.splitlines():
        if re.match(rf"^  {section}:\s*$", line):
            inside = True
        elif inside and re.match(r"^\s{4}- ", line):
            entries.append(line.split("#")[0].strip()[2:])
        elif inside and line.strip():
            break
    logger.debug("meta.yaml %s requirements: %s", section, entries)
    return dict(split(entry) for entry in entries)


def main():
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    recipe = (ROOT / ".github/conda/meta.yaml").read_text(encoding="utf-8")
    run = recipe_section(recipe, "run")
    host = recipe_section(recipe, "host")

    expected = {}
    for requirement in pyproject["project"]["dependencies"]:
        name, specifier = split(requirement)
        expected[CONDA_NAMES.get(name, name)] = specifier

    problems = []
    for name, specifier in expected.items():
        if name not in run:
            problems.append(f"{name}{specifier} is in pyproject.toml but not meta.yaml run")
        elif run[name] != specifier:
            problems.append(f"{name}: pyproject.toml has '{specifier}', meta.yaml run has '{run[name]}'")
    for name in run:
        if name not in expected and name != "python" and name not in host:
            problems.append(f"{name} is in meta.yaml run but not pyproject.toml or host")

    for problem in problems:
        logger.error(problem)
    if problems:
        return 1
    logger.info("meta.yaml run requirements match pyproject.toml (%d packages)", len(expected))
    return 0


if __name__ == "__main__":
    sys.exit(main())
