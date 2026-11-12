"""Where every file of the R&D baseline goes.

The baseline (the commit ``BASELINE`` below) kept 165 design notes in ``docs/``,
1,174 verification records in ``docs/verification/`` and 861 experiment files in
``tools/design-experiments/``, each in one flat directory. This module maps each
baseline path to its place in the organized tree. It holds data and pure
functions only; ``reorganize.py`` applies it.

Target layout::

    docs/requirements/   what the product must do (user decisions, acceptance)
    docs/specs/<area>/   adopted specifications, the canon for implementation
    docs/plans/          the start review: packages, order, risks
    rd/planning/         plans and ledgers of the pre-implementation work (history)
    docs/engineering/    rules and lessons
    docs/papers/         technical report
    rd/<domain>/<family>/            design notes of one line of investigation
    rd/<domain>/<family>/records/    its verification records
    rd/<domain>/<family>/experiments/  its scripts
    rd/<domain>/<family>/fixtures/   its input decks
"""

from __future__ import annotations

import re
from pathlib import PurePosixPath

# The commit that holds the R&D state before reorganization. It is named by its
# hash, not by a tag: the commit is in the history of main, and a hash cannot be moved.
BASELINE = "51bd5b544e7b821cdbbbb01ad8be4f53ca0e6a5b"
# Written when the baseline was redacted; it travels with the baseline.
REDACTION_RECORD = "redaction-record.json"
OLD_RECORDS = "docs/verification"
OLD_EXPERIMENTS = "tools/design-experiments"

# Leading verbs of script names; the topic follows them.
VERBS = (
    "check|audit|run|prepare|assess|close|diagnose|report|probe|sync|summarize|analyze|launch|"
    "export|finalize|fetch|certify|compare|build|save|plot|extract|inspect|reassess|reevaluate|"
    "publish|update"
)

# (family, pattern). The first pattern that matches the topic wins, so the
# specific rows stand above the general ones.
FAMILY_RULES: list[tuple[str, str]] = [
    # --- coupled problem, time integration, stability
    (
        "coupled/j1",
        (
            r"^(j1|coupled|age|harmonic-age|material-birth|so3|static-branch|flutter-search|"
            r"pg-uvlm|uvlm-age|load-release-initialization)"
        ),
    ),
    # --- acceptance and start review
    ("acceptance/s10", r"^s10"),
    ("acceptance/plan", r"^(practical-acceptance|acceptance-case|implementation-readiness)"),
    # --- build and numerical foundation
    ("foundation/build", r"^(build-2026|cpp23|2026-10-29)"),
    ("foundation/numerics", r"^(numerical|io-dependencies|linear-solve)"),
    # --- structure: main shell candidates, in the order they were studied
    ("structure/main-shell/ancf-3423", r"^(3423|ancf-3423|chrono-shell|eas-layer)"),
    (
        "structure/main-shell/q9-dynamics",
        (
            r"^(q9-(event|integrator|load|material|multielement|nonlinear|prestress|transient)|"
            r"prestress-equilibrium|nonlinear-step|transient|multielement-prestress)"
        ),
    ),
    ("structure/main-shell/q9-curved", r"^q9-(curved|developable|eas21)"),
    ("structure/main-shell/q9-edge-trace", r"^q9-(edge|trace)"),
    (
        "structure/main-shell/q9-plate",
        r"^(3833|ancf-3833|q9|plate-modal|modal-reference|ancf-shell-candidate|curved-shell-comparison)",
    ),
    ("structure/main-shell/q3", r"^(q3|developable-general|general-bending)"),
    ("structure/main-shell/mq3", r"^(mq3|curved-shell-(moment|joint))"),
    ("structure/main-shell/gsp-q3", r"^(gsp-q3|main-shell-formulation)"),
    ("structure/main-shell/mixed-order-150", r"^mixed-(150|order)"),
    (
        "structure/main-shell/q4-projection",
        (
            r"^(qem|mid-projection|midsurface-projection|covariant|geometric-(projection|mid)|tensor|"
            r"tying|q4-continuity|main-shell-projection)"
        ),
    ),
    ("structure/main-shell/mortar", r"^(mortar|joint-mortar)"),
    ("structure/main-shell/hybrid-gradient", r"^hybrid-gradient"),
    ("structure/main-shell/polynomial-trace", r"^polynomial-(trace|reference)"),
    (
        "structure/main-shell/reference-work",
        r"^(reference-work|curved-work|(stable-)?curvature-(profiles|work))",
    ),
    ("structure/main-shell/moment-enriched", r"^moment-(enriched|position|recovery|exact)"),
    ("structure/main-shell/mixed-mean-stress", r"^mixed-mean"),
    (
        "structure/main-shell/geometry-compatible",
        r"^(geometry-compatible|c1-work|gc-|coefficient-compatible)",
    ),
    (
        "structure/main-shell/paired-geometry",
        r"^(paired|pg-|mass-paired|coefficient-paired|exact-paired|general-paired|swept-taper)",
    ),
    (
        "structure/main-shell/review",
        r"^(main-shell|structural-progress-control|structural-model-adoption)",
    ),
    # --- structure: the rest
    ("structure/material", r"^(material|initial-strain)"),
    ("structure/beam", r"^beam"),
    ("structure/connection", r"^connection"),
    ("structure/integration", r"^structural"),
    ("structure/conventions", r"^(thickness-coordinate|reference-geometry)"),
    # --- aerodynamics
    (
        "aero/avl-reference",
        (
            r"^(avl|fixed-band|symmetric-reference|spanwise-load|steady-horseshoe|field-representation|"
            r"reference-avl)"
        ),
    ),
    (
        "aero/finite-wing",
        r"^(finite-wing|trefftz|drag|induced-drag|steady-drag|edge-load|uvlm-drag|uvlm-steady)",
    ),
    ("aero/multi-surface", r"^(multiple-surface|uvlm-multiple|uvlm-global)"),
    (
        "aero/unsteady-reference",
        (
            r"^(unsteady|frequency|high-frequency|ideal|pitch-frequency|symmetric-heave|ptera|"
            r"aerosandbox|aerodynamic-reference|analytic-reference|reference-campaign|"
            r"uvlm-(harmonic|heave|pitch|frequency|adaptive|error|independent|reference|analytic))"
        ),
    ),
    (
        "aero/kernel",
        r"^(uvlm-(induction|near|rm|dual13|edge|surface|geometry)|vortex-core|aerogrid)",
    ),
    ("aero/adoption", r"^uvlm-(adoption|integrated)"),
    ("aero/wake-loads", r"^uvlm"),
]

# Design notes that are not R&D history: they move into docs/, by exact name.
DOCS_PLACEMENT: dict[str, str] = {
    # requirements and acceptance: user decisions and what "done" means
    "product-requirements.md": "docs/requirements",
    "project-start.md": "docs/requirements",
    "acceptance-criteria.md": "docs/requirements",
    "performance-plan.md": "docs/requirements",
    "practical-acceptance-plan.md": "docs/requirements",
    "requirements-evidence-map.md": "rd/planning",
    "verification-resources.md": "docs/requirements",
    # adopted specifications
    "coupled-problem-definition.md": "docs/specs",
    "numerical-conventions.md": "docs/specs",
    "structural-adopted-specification.md": "docs/specs/structure",
    "main-shell-adopted-specification.md": "docs/specs/structure",
    "material-shell-adopted-specification.md": "docs/specs/structure",
    "beam-adopted-specification.md": "docs/specs/structure",
    "connection-adopted-specification.md": "docs/specs/structure",
    "uvlm-integrated-specification.md": "docs/specs/aero",
    "uvlm-induction-adopted-specification.md": "docs/specs/aero",
    "uvlm-wake-adopted-specification.md": "docs/specs/aero",
    "coupled-adopted-specification.md": "docs/specs/coupled",
    "j1-initialization-output-derivative-addendum.md": "docs/specs/coupled",
    "j1-discrete-equilibrium-and-chart-contract.md": "docs/specs/coupled",
    "j1-static-continuation-error-contract.md": "docs/specs/coupled",
    "coupled-time-stability-specification.md": "docs/specs/coupled",
    "coupled-initialization-residual-contract.md": "docs/specs/coupled",
    "coupled-state-search-contract.md": "docs/specs/coupled",
    "numerical-foundation-adoption-specification.md": "docs/specs/foundation",
    "io-execution-adoption-specification.md": "docs/specs/foundation",
    "linear-solve-accuracy-contract.md": "docs/specs/foundation",
    # Defines MAT-D04..07, which the adopted material and shell specifications cite.
    "material-layer-adoption-contract.md": "docs/specs/structure",
    # Defines MAT-D01..03 (reference axes, consistent tangent, one material-point evaluation).
    "material-initial-strain-contract.md": "docs/specs/structure",
    # the start review is the current plan; the other plans and ledgers are history
    # The page that tracked the pre-implementation work. docs/roadmap.md is a new page.
    "roadmap.md": "rd/planning/preimplementation-history.md",
    "preimplementation-plan.md": "rd/planning",
    "implementation-start-review.md": "docs/plans",
    "implementation-backlog.md": "rd/planning",
    # The ledger defines decision IDs (AE-01..07 and others) that the specifications cite.
    "design-decisions.md": "docs/specs",
    "design-contracts-and-invariants.md": "rd/planning",
    "structural-model-adoption-review.md": "rd/planning",
    "structural-integration-after-shell.md": "rd/planning",
    "restart-and-distribution.md": "rd/planning",
    "uvlm-ancf-oss-development-plan.md": "rd/planning",
    "uvlm-ancf-software-design.md": "rd/planning",
    "build.md": "docs",
}

# Whole baseline directories that move as they are.
TREE_MOVES: list[tuple[str, str]] = [
    ("docs/engineering/knowledge-validation.md", "rd/planning/knowledge-validation.md"),
    ("docs/engineering", "docs/engineering"),
    ("docs/papers", "docs/papers"),
    ("output/pdf", "docs/papers/pdf"),
    # The scripts that generated the technical report; frozen with the report they made.
    ("tools/papers", "rd/structure/main-shell/paired-geometry/report-tools"),
    ("tools/toolchains", "tools/toolchains"),
    ("tools/skills", ".agents/skills"),
    (OLD_EXPERIMENTS + "/cpp-numerics", "rd/foundation/numerics/experiments/cpp-numerics"),
    (OLD_EXPERIMENTS + "/cpp-io-binding", "rd/foundation/numerics/experiments/cpp-io-binding"),
    (OLD_RECORDS + "/fixtures", "rd/aero/avl-reference/fixtures"),
    (OLD_RECORDS + "/figures", "rd/aero/avl-reference/figures"),
    (OLD_EXPERIMENTS + "/requirements-reference.txt", "rd/requirements-reference.txt"),
]

# Literal edits to single scripts where a path test cannot be rewritten by rule.
SCRIPT_PATCHES: dict[str, list[tuple[str, str]]] = {
    OLD_EXPERIMENTS + "/build_acceptance_case_catalog.py": [
        ("k.startswith('tools/design-experiments/')", "'/experiments/' in k"),
    ],
}

# (record, key) pairs whose hash names a copy kept under build/, not the
# repository file with the same baseline content. Found by running the audit
# scripts in the baseline and in the organized tree.
HASH_KEEP: set[tuple[str, str]] = {
    (OLD_RECORDS + "/avl-field-overlay-decimal-result.json", "avl_field_overlay.py"),
}

# Pages rewritten by hand after the reorganization. They are outside the claim
# that the tree equals the transform of the baseline.
MAINTAINED: set[str] = {
    "README.md",
    "AGENTS.md",
    "docs/build.md",
    "docs/engineering/README.md",
    # A new page at the path the pre-implementation roadmap had; the old page
    # is rd/planning/preimplementation-history.md.
    "docs/roadmap.md",
}

# Files of the baseline that were removed afterwards, with the reason. They stay
# in the baseline tag. Links to them from the records are left as they are.
REMOVED: dict[str, str] = {
    "tools/toolchains/README.md": (
        "Isolated GCC 14 toolchain of the 2026-10-29 Linux check; "
        "not used by the product build or CI"
    ),
    "tools/toolchains/linux-gcc14.explicit.txt": "same",
    "tools/toolchains/linux-gcc14.packages.json": "same",
}

# (record, cited path) pairs whose file is, by its nature, neither in the
# repository nor in the evidence index, with the reason. check_tree.py counts
# these and fails on every other cited path that it cannot find.
_AVL_LINK_INPUT = (
    "Static library that the AVL build makes from the AVL source; "
    "an intermediate of the link, not kept after the executables were built"
)
_CHRONO_SOURCE = (
    "File of Project Chrono at the revision the record names, with its URL; "
    "the path is relative to that repository"
)
HASH_EXTERNAL: dict[tuple[str, str], str] = {
    **{
        (f"{OLD_RECORDS}/{record}.json", f"build/reference-avl/source-baseline/{library}"): (
            _AVL_LINK_INPUT
        )
        for record in (
            "avl-span1280-capacity-ks-preparation",
            "avl-span1280-capacity-preparation",
            "avl-symmetric-solve-preparation",
        )
        for library in ("plotlib/libPlt_gDP.a", "eispack/libeispack.a")
    },
    **{
        (f"{OLD_RECORDS}/chrono-shell-source-manifest.json", name): _CHRONO_SOURCE
        for name in (
            "LICENSE",
            "src/chrono/fea/ChElementShellANCF_3423.cpp",
            "src/chrono/fea/ChElementShellANCF_3423.h",
            "src/chrono/fea/ChElementShellANCF_3833.cpp",
            "src/chrono/fea/ChElementShellANCF_3833.h",
            "src/chrono/fea/ChNodeFEAxyzDD.h",
        )
    },
}

# Files outside the R&D material keep their place.
UNMOVED_PREFIXES = (
    ".gitattributes",
    ".gitignore",
    "AGENTS.md",
    "LICENSE",
    "README.md",
    "cmake/",
    "include/",
    "src/",
    "tests/",
)

_RULES = [(family, re.compile(pattern)) for family, pattern in FAMILY_RULES]
_VERB = re.compile(rf"^(?:{VERBS})-")


def topic(name: str) -> str:
    """File name -> the words that name its line of investigation."""
    stem = name.rsplit(".", 1)[0].replace("_", "-").lower()
    previous = None
    while previous != stem:  # e.g. run_check_x
        previous = stem
        stem = _VERB.sub("", stem)
    return stem


def family_of(name: str) -> str | None:
    key = topic(name)
    for family, pattern in _RULES:
        if pattern.search(key):
            return family
    return None


def new_path(old: str) -> str | None:
    """Baseline path -> organized path. None means no rule covers the file."""
    if old == REDACTION_RECORD:
        return "rd/" + REDACTION_RECORD
    if old == ".gitattributes" or old.startswith(UNMOVED_PREFIXES):
        return old
    for source, target in TREE_MOVES:
        if old == source or old.startswith(source + "/"):
            return target + old[len(source) :]
    path = PurePosixPath(old)
    parent, name = path.parent.as_posix(), path.name
    if parent == "docs":
        if name in DOCS_PLACEMENT:
            place = DOCS_PLACEMENT[name]
            return place if place.endswith(".md") else f"{place}/{name}"
        family = family_of(name)
        return f"rd/{family}/{name}" if family else None
    if parent == OLD_RECORDS:
        family = family_of(name)
        return f"rd/{family}/records/{name}" if family else None
    if parent == OLD_EXPERIMENTS:
        family = family_of(name)
        return f"rd/{family}/experiments/{name}" if family else None
    return None
