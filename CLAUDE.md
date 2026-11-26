# CLAUDE.md

Instructions for Claude Code in this repository. They add to `AGENTS.md`,
which stays the shared rulebook for all coding agents:

@AGENTS.md

Where the two differ, this file wins for Claude. In particular, the Codex
model names and the `.codex/agents/` role files in `AGENTS.md` and `README.md`
do not apply to Claude. Use the model routing below instead.

## Project in one paragraph

CMS dijet analysis of pPb collisions at sqrt(s_NN) = 8.16 TeV, with pp
references. `processForestSimple.C` (about 4.6k lines) is both a ROOT macro and
the source of the compiled executable `build/processForestSimple`. It reads
HiForest trees through `TChain`s, applies event selection, JEC/JER and jet ID,
and writes histograms to one ROOT file. `hist_analysis/` (PyROOT, RooUnfold,
Jupyter, a C++17 regularization scan in `hist_analysis/cpp/`) turns those
histograms into closures, systematics, unfolded results and comparisons with
theory. `processing/` holds the HTCondor workflow. The measurement plan and its
approval-gated steps live in `steps/README.md`. `README.md` is the source of
truth for commands, arguments and validation. Do not restate it here; read it.

## Do not invent facts

- Every statement about physics, code behavior, data, or results must be
  backed by something you read or ran: a file and line, a command output, a
  ROOT key, a cited paper or official CMS documentation. Cite it.
- If you cannot verify something, say "not verified" and say what would verify
  it. Never fill gaps with plausible numbers, branch names, histogram keys,
  correction-file names, cross sections, luminosities, or paper results.
- Never claim a build, test, or notebook passed unless you ran it in this
  session. Report skipped checks and the exact reason (for example, input
  forests or ROOT not available).
- Undocumented analysis conventions are the user's decision. Ask one focused
  question instead of guessing (this repeats `AGENTS.md`, on purpose).

## Work in agentic mode

Work autonomously through plan, critique, implement, validate, report. Do not
stop after each small step to ask permission for reversible work, but stop
and ask before anything listed under "Always ask first".

1. **Understand.** Read the relevant code, `README.md` sections, the matching
   `steps/` document, and `git status --short`.
2. **Plan.** Write a short plan. For a physics-relevant change, include what
   `AGENTS.md` requires: observable, provenance, units and frame, selections
   and weights, binning and flow, normalization, uncertainty prescription,
   expected changes, acceptance criteria.
3. **Challenge the plan (mandatory, every time).** Before implementing, launch
   an independent subagent to critically assess the plan. See the next section.
4. **Implement** the smallest change that satisfies the plan.
5. **Validate** with the smallest applicable checks from `README.md`. For
   scientifically significant changes, numerical validation is done by a
   separate agent, not by the implementer.
6. **Report** with the handoff in `AGENTS.md`: files touched, behavior change,
   checks run or skipped, remaining uncertainty, expected physics or output
   differences.

### Independent critical assessment of plans

Always run a separate, fresh subagent (not a fork of the current context) to
challenge any proposed plan before it is executed.

- Give the critic the user's request, the relevant source files, and the
  requirements first, and your plan second, so it forms its own expectations
  before reading your reasoning.
- Instruct it to be adversarial: find wrong assumptions, physics errors
  (frames, eta shift and beam direction, weights, normalization, flow bins,
  error propagation, miss/fake definitions), unprotected changes to the
  guardrails in `AGENTS.md`, missing validation, and simpler alternatives.
- The critic is read-only. It must not edit files.
- It must cite file and line evidence for each objection, and say plainly when
  it finds no material problem.
- Resolve each objection by evidence, or revise the plan. Objections that turn
  on an analysis choice go to the user. Summarize the critique and how it was
  resolved in your reply.
- Scale the critic's model to the risk (see routing). A documentation plan
  can be checked by a small model; a change to selections, corrections,
  responses or unfolding gets the strongest model.

## Model routing

Use the cheapest model that can do the task well. Pass the model explicitly
when launching a subagent.

| Task | Model |
|---|---|
| Documentation writing, README or comment edits, summarizing logs, listing files or keys, simple searches | `haiku` |
| Routine code edits, tests, notebook plumbing, batch-script changes, plan critique for low-risk work | `sonnet` |
| Physics-relevant changes (selections, JEC/JER, weights, matching, response, unfolding, systematics), physics review, numerical validation design, hard debugging | `opus` |

Keep at most three subagents running at once, give each its own files
(one writer per file), and do not let the implementer certify its own work.

## Tools

Use the tools that senior C++/Python developers and HEP physicists use, and
prefer the ones this repository already relies on. Check that a tool is
installed (`which`, `root-config --version`, `py-env/bin/python -m pip list`)
before relying on it, and do not assume it exists.

Already used in this repository:

- **ROOT / PyROOT** (macros, `TChain`, histograms, `rootcling` dictionaries),
  **RooUnfold** (Bayesian unfolding), **mplhep** and the shared style in
  `hist_analysis/python/root_style.py`.
- **CMake + CTest** (`cmake -S . -B build`, `ctest --test-dir build`), with
  Release first and Debug to reproduce failures. CI uses
  `rootproject/root:6.26.10-ubuntu22.04` (`.github/workflows/ci.yml`).
- **Python `unittest`** under `tests/`, `py_compile`/`compileall` syntax checks,
  `bash -n` for shell scripts, Jupyter run from the repository root with
  `py-env/bin/python`.
- **git** for inspection (`git status`, `git diff`, `git log -p`, `git blame`).
- **HTCondor**, EOS and XRootD for batch work, always with `--dry-run` first.
- **ROOT MCP server** (`.codex/run-root-mcp.sh`, see `README.md`) for ROOT-file
  inspection, if it is configured for Claude.

Standard tools to propose when they fit the problem: `gdb`/`lldb` for crashes,
AddressSanitizer/UBSan or Valgrind for memory errors, `clang-tidy` for static
analysis, `perf`/Instruments or ROOT's `TStopwatch` for performance, and
`rootls`/`rootprint` or `uproot` for quick file inspection.

**Ask the user before introducing any tool, dependency, or formatter that the
repository does not already use**, and say why it helps. Do not run a
formatter over `processForestSimple.C` (formatting-only diffs hide physics
changes, see `AGENTS.md`).

### Context7 MCP for library documentation

Use the Context7 MCP server to fetch current documentation before writing or
changing code that calls an external library API (ROOT, RooUnfold, mplhep,
matplotlib, NumPy, uproot, CMake, HTCondor Python bindings, pytest, and so on).
Resolve the library ID first, then query the docs for the specific class or
function. Prefer this over memory, especially for ROOT, where APIs differ
between versions; check the local ROOT version with `root-config --version`.

If Context7 is not available in the session, say so, fall back to the
official documentation (root.cern, the library's own docs) or the installed
headers and source, and do not guess signatures. To add it to Claude Code:
`claude mcp add --transport http context7 https://mcp.context7.com/mcp`.

## Always ask first

These need the user's explicit approval for that specific operation:

- any change to the protected physics behavior listed in `AGENTS.md`
  (selections, corrections, weights, binning, histogram keys, ROOT output
  structure, response/fake/miss definitions);
- HTCondor submission, job cancellation, deleting or replacing validated
  outputs, merges that drop inputs;
- new dependencies or tools, and edits to `aux_files/` (treat it as data);
- `git commit`, `push`, rebases, or anything that rewrites history.

## Review reports

Save review findings under `review/`, as `AGENTS.md` requires. Name files that
Claude writes with the prefix `Claude_` (for example
`review/Claude_processForestSimple_review.md`). Reports in `review/` are local
artifacts and are not committed (`review/` is in `.gitignore`).
