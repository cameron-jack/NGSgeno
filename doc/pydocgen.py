#!/usr/bin/env python3
"""
pydocgen.py — Python documentation generator with call graph analysis.

Usage:
    python pydocgen.py <source_file_or_dir> [--output docs.html] [--title "My Project"]

Extracts:
  - Module, class, function, and global docstrings
  - Caller/callee relationships (static analysis via AST)
  - Type annotations
  - Produces a self-contained HTML documentation site
"""

import ast
import os
import sys
import json
import argparse
import textwrap
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Optional


# ─── Data models ─────────────────────────────────────────────────────────────

@dataclass
class ParamInfo:
    name: str
    annotation: Optional[str] = None
    default: Optional[str] = None


@dataclass
class FunctionDoc:
    name: str
    qualified_name: str          # e.g. MyClass.my_method
    module: str
    kind: str                    # "function" | "method" | "classmethod" | "staticmethod"
    docstring: Optional[str]
    params: list[ParamInfo]
    returns: Optional[str]
    lineno: int
    calls: list[str] = field(default_factory=list)      # names this function calls
    callers: list[str] = field(default_factory=list)    # filled in post-pass


@dataclass
class ClassDoc:
    name: str
    qualified_name: str
    module: str
    docstring: Optional[str]
    bases: list[str]
    methods: list[FunctionDoc]
    lineno: int


@dataclass
class GlobalDoc:
    name: str
    qualified_name: str
    module: str
    annotation: Optional[str]
    value: Optional[str]
    docstring: Optional[str]      # comment/string immediately after assignment
    lineno: int


@dataclass
class ModuleDoc:
    name: str
    path: str
    docstring: Optional[str]
    functions: list[FunctionDoc]
    classes: list[ClassDoc]
    globals: list[GlobalDoc]


# ─── AST helpers ─────────────────────────────────────────────────────────────

def annotation_to_str(node) -> Optional[str]:
    if node is None:
        return None
    return ast.unparse(node)


def default_to_str(node) -> Optional[str]:
    if node is None:
        return None
    return ast.unparse(node)


def collect_calls(node: ast.AST) -> list[str]:
    """Return all function/method names called within an AST subtree."""
    calls = []
    for child in ast.walk(node):
        if isinstance(child, ast.Call):
            func = child.func
            if isinstance(func, ast.Name):
                calls.append(func.id)
            elif isinstance(func, ast.Attribute):
                calls.append(func.attr)
    return list(dict.fromkeys(calls))   # deduplicate, preserve order


def parse_params(args: ast.arguments) -> list[ParamInfo]:
    params = []
    # Align defaults: they apply to the last N args
    all_args = args.posonlyargs + args.args
    defaults_offset = len(all_args) - len(args.defaults)

    for i, arg in enumerate(all_args):
        default_idx = i - defaults_offset
        default = default_to_str(args.defaults[default_idx]) if default_idx >= 0 else None
        params.append(ParamInfo(
            name=arg.arg,
            annotation=annotation_to_str(arg.annotation),
            default=default,
        ))

    if args.vararg:
        params.append(ParamInfo(name=f"*{args.vararg.arg}",
                                annotation=annotation_to_str(args.vararg.annotation)))
    for kw in args.kwonlyargs:
        params.append(ParamInfo(name=kw.arg, annotation=annotation_to_str(kw.annotation)))
    if args.kwarg:
        params.append(ParamInfo(name=f"**{args.kwarg.arg}",
                                annotation=annotation_to_str(args.kwarg.annotation)))
    return params


def kind_of_method(node: ast.FunctionDef) -> str:
    for dec in node.decorator_list:
        if isinstance(dec, ast.Name):
            if dec.id == "classmethod":
                return "classmethod"
            if dec.id == "staticmethod":
                return "staticmethod"
        elif isinstance(dec, ast.Attribute) and dec.attr in ("classmethod", "staticmethod"):
            return dec.attr
    return "method"


# ─── Module parser ───────────────────────────────────────────────────────────

def parse_module(path: Path, base_dir: Path) -> ModuleDoc:
    source = path.read_text(encoding="utf-8", errors="replace")
    try:
        tree = ast.parse(source, filename=str(path))
    except SyntaxError as e:
        print(f"  Warning: syntax error in {path}: {e}", file=sys.stderr)
        return ModuleDoc(name=path.stem, path=str(path),
                         docstring=None, functions=[], classes=[], globals=[])

    rel = path.relative_to(base_dir)
    module_name = ".".join(rel.with_suffix("").parts)

    module_doc = ast.get_docstring(tree)
    functions: list[FunctionDoc] = []
    classes: list[ClassDoc] = []
    globals_: list[GlobalDoc] = []

    # Track which names are classes/functions (for global detection)
    defined_names = set()

    for node in ast.iter_child_nodes(tree):
        # ── Top-level functions ───────────────────────────────────────────
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            defined_names.add(node.name)
            fdoc = FunctionDoc(
                name=node.name,
                qualified_name=f"{module_name}.{node.name}",
                module=module_name,
                kind="function",
                docstring=ast.get_docstring(node),
                params=parse_params(node.args),
                returns=annotation_to_str(node.returns),
                lineno=node.lineno,
                calls=collect_calls(node),
            )
            functions.append(fdoc)

        # ── Classes ───────────────────────────────────────────────────────
        elif isinstance(node, ast.ClassDef):
            defined_names.add(node.name)
            bases = [ast.unparse(b) for b in node.bases]
            methods: list[FunctionDoc] = []

            for item in ast.iter_child_nodes(node):
                if isinstance(item, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    mdoc = FunctionDoc(
                        name=item.name,
                        qualified_name=f"{module_name}.{node.name}.{item.name}",
                        module=module_name,
                        kind=kind_of_method(item),
                        docstring=ast.get_docstring(item),
                        params=parse_params(item.args),
                        returns=annotation_to_str(item.returns),
                        lineno=item.lineno,
                        calls=collect_calls(item),
                    )
                    methods.append(mdoc)

            classes.append(ClassDoc(
                name=node.name,
                qualified_name=f"{module_name}.{node.name}",
                module=module_name,
                docstring=ast.get_docstring(node),
                bases=bases,
                methods=methods,
                lineno=node.lineno,
            ))

        # ── Module-level assignments (globals) ────────────────────────────
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
                name = node.target.id
                annotation = annotation_to_str(node.annotation)
                value = ast.unparse(node.value) if node.value else None
            elif isinstance(node, ast.Assign):
                # Only handle simple single-name targets
                if len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
                    name = node.targets[0].id
                    annotation = None
                    value = ast.unparse(node.value)
                else:
                    continue
            else:
                continue

            if name.startswith("_"):
                continue   # skip private/dunder globals

            defined_names.add(name)
            globals_.append(GlobalDoc(
                name=name,
                qualified_name=f"{module_name}.{name}",
                module=module_name,
                annotation=annotation,
                value=value[:80] + "…" if value and len(value) > 80 else value,
                docstring=None,
                lineno=node.lineno,
            ))

    return ModuleDoc(
        name=module_name,
        path=str(path),
        docstring=module_doc,
        functions=functions,
        classes=classes,
        globals=globals_,
    )


# ─── Project-level analysis ───────────────────────────────────────────────────

def collect_python_files(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    files = []
    for p in sorted(path.rglob("*.py")):
        # skip hidden dirs, __pycache__, venv, etc.
        parts = p.parts
        if any(part.startswith(".") or part in ("__pycache__", "venv", ".venv", "node_modules")
               for part in parts):
            continue
        files.append(p)
    return files


def resolve_call_graph(modules: list[ModuleDoc]) -> dict[str, list[str]]:
    """
    Build a reverse mapping: callee_short_name → list of caller qualified_names.
    Returns callers dict keyed by short function name.
    """
    # All known qualified names and their short names
    all_funcs: dict[str, FunctionDoc] = {}
    for m in modules:
        for f in m.functions:
            all_funcs[f.qualified_name] = f
            all_funcs[f.name] = f          # also index by short name
        for c in m.classes:
            for mth in c.methods:
                all_funcs[mth.qualified_name] = mth
                all_funcs[mth.name] = mth

    callers_map: dict[str, list[str]] = defaultdict(list)

    for m in modules:
        for caller in m.functions:
            for callee_name in caller.calls:
                callers_map[callee_name].append(caller.qualified_name)
        for c in m.classes:
            for caller in c.methods:
                for callee_name in caller.calls:
                    callers_map[callee_name].append(caller.qualified_name)

    # Attach callers back onto each function
    for m in modules:
        for f in m.functions:
            seen = callers_map.get(f.name, []) + callers_map.get(f.qualified_name, [])
            f.callers = list(dict.fromkeys(seen))
        for c in m.classes:
            for mth in c.methods:
                seen = callers_map.get(mth.name, []) + callers_map.get(mth.qualified_name, [])
                mth.callers = list(dict.fromkeys(seen))

    return callers_map


# ─── JSON serialisation for the HTML template ─────────────────────────────────

def to_json(modules: list[ModuleDoc]) -> str:
    def ser(obj):
        if hasattr(obj, "__dataclass_fields__"):
            return obj.__dict__
        raise TypeError(repr(obj))
    return json.dumps([m.__dict__ for m in modules], default=ser)


# ─── HTML generation ──────────────────────────────────────────────────────────

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>{{TITLE}}</title>
<style>
  @import url('https://fonts.googleapis.com/css2?family=DM+Mono:ital,wght@0,300;0,400;0,500;1,400&family=Fraunces:ital,opsz,wght@0,9..144,300;0,9..144,600;1,9..144,300&display=swap');

  :root {
    --bg: #0e0f11;
    --surface: #16181d;
    --surface2: #1e2028;
    --border: #2a2d38;
    --accent: #e8c84a;
    --accent2: #5b8af0;
    --accent3: #8b5cf6;
    --text: #d4d8e8;
    --text-dim: #7b8290;
    --text-bright: #f0f2ff;
    --green: #34d399;
    --red: #f87171;
    --orange: #fb923c;
    --sidebar-w: 280px;
  }

  * { box-sizing: border-box; margin: 0; padding: 0; }

  body {
    font-family: 'DM Mono', monospace;
    background: var(--bg);
    color: var(--text);
    display: flex;
    min-height: 100vh;
    font-size: 14px;
    line-height: 1.6;
  }

  /* ── Sidebar ── */
  #sidebar {
    width: var(--sidebar-w);
    min-width: var(--sidebar-w);
    background: var(--surface);
    border-right: 1px solid var(--border);
    height: 100vh;
    position: sticky;
    top: 0;
    overflow-y: auto;
    display: flex;
    flex-direction: column;
  }
  #sidebar-header {
    padding: 20px 18px 14px;
    border-bottom: 1px solid var(--border);
  }
  #sidebar-header h1 {
    font-family: 'Fraunces', serif;
    font-size: 17px;
    color: var(--accent);
    font-weight: 600;
    letter-spacing: .3px;
  }
  #sidebar-header p {
    font-size: 12px;
    color: var(--text-dim);
    margin-top: 3px;
  }
  #search-box {
    margin: 10px 12px;
    padding: 7px 10px;
    width: calc(100% - 24px);
    background: var(--surface2);
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--text);
    font-family: 'DM Mono', monospace;
    font-size: 12px;
    outline: none;
  }
  #search-box:focus { border-color: var(--accent2); }
  #search-box::placeholder { color: var(--text-dim); }

  .nav-module { margin-bottom: 2px; }
  .nav-module-header {
    padding: 7px 12px;
    font-size: 12px;
    color: var(--accent);
    letter-spacing: .5px;
    text-transform: uppercase;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 6px;
    user-select: none;
  }
  .nav-module-header:hover { background: var(--surface2); }
  .nav-module-header .chevron { transition: transform .2s; font-size: 14px; }
  .nav-module-header.collapsed .chevron { transform: rotate(-90deg); }
  .nav-items { padding-left: 4px; }
  .nav-item {
    padding: 4px 12px 4px 18px;
    font-size: 12px;
    color: var(--text-dim);
    cursor: pointer;
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    display: flex;
    align-items: center;
    gap: 5px;
    border-radius: 4px;
    margin: 1px 6px;
    text-decoration: none;
  }
  .nav-item:hover { background: var(--surface2); color: var(--text); }
  .nav-item.active { background: rgba(232,200,74,.12); color: var(--accent); }
  .nav-item .kind-badge { font-size: 9px; opacity: .6; }

  /* ── Main content ── */
  #main {
    flex: 1;
    padding: 40px 52px;
    max-width: 960px;
  }

  .page-title {
    font-family: 'Fraunces', serif;
    font-size: 38px;
    color: var(--text-bright);
    font-weight: 300;
    margin-bottom: 6px;
    letter-spacing: -.5px;
  }
  .page-subtitle {
    color: var(--text-dim);
    font-size: 13px;
    margin-bottom: 48px;
  }

  /* ── Module section ── */
  .module-section { margin-bottom: 64px; }
  .module-name {
    font-family: 'Fraunces', serif;
    font-size: 24px;
    color: var(--accent);
    font-weight: 600;
    padding-bottom: 8px;
    border-bottom: 1px solid var(--border);
    margin-bottom: 18px;
    display: flex;
    align-items: baseline;
    gap: 10px;
  }
  .module-path { font-size: 12px; color: var(--text-dim); font-family: 'DM Mono', monospace; font-weight:300; }
  .module-docstring {
    background: var(--surface);
    border-left: 3px solid var(--accent2);
    padding: 12px 16px;
    border-radius: 0 6px 6px 0;
    margin-bottom: 28px;
    color: var(--text);
    font-size: 13px;
    white-space: pre-wrap;
  }

  /* ── Entry cards ── */
  .entry-card {
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 10px;
    margin-bottom: 18px;
    overflow: hidden;
    transition: border-color .15s;
  }
  .entry-card:hover { border-color: #3a3f55; }
  .entry-card:target { border-color: var(--accent2); }

  .entry-header {
    padding: 14px 18px;
    cursor: pointer;
    display: flex;
    align-items: flex-start;
    gap: 12px;
    user-select: none;
  }
  .entry-header:hover { background: rgba(255,255,255,.02); }

  .entry-kind {
    font-size: 10px;
    padding: 2px 7px;
    border-radius: 4px;
    font-weight: 500;
    letter-spacing: .4px;
    text-transform: uppercase;
    margin-top: 2px;
    flex-shrink: 0;
  }
  .kind-function  { background: rgba(91,138,240,.18); color: var(--accent2); }
  .kind-method    { background: rgba(139,92,246,.18); color: var(--accent3); }
  .kind-class     { background: rgba(52,211,153,.15); color: var(--green); }
  .kind-classmethod { background: rgba(251,146,60,.15); color: var(--orange); }
  .kind-staticmethod { background: rgba(248,113,113,.15); color: var(--red); }
  .kind-global    { background: rgba(232,200,74,.12); color: var(--accent); }

  .entry-sig {
    flex: 1;
    min-width: 0;
  }
  .entry-name {
    font-family: 'DM Mono', monospace;
    font-size: 15px;
    color: var(--text-bright);
    font-weight: 500;
  }
  .entry-return {
    font-size: 12px;
    color: var(--text-dim);
    margin-top: 2px;
  }
  .entry-return span { color: var(--accent3); }

  .entry-body {
    border-top: 1px solid var(--border);
    padding: 16px 18px;
    display: none;
  }
  .entry-body.open { display: block; }

  .docstring-text {
    font-size: 13px;
    color: var(--text);
    white-space: pre-wrap;
    margin-bottom: 14px;
    line-height: 1.7;
  }
  .no-doc { color: var(--text-dim); font-style: italic; font-size: 12px; margin-bottom: 14px; }

  /* Params table */
  .params-table {
    width: 100%;
    border-collapse: collapse;
    font-size: 12px;
    margin-bottom: 14px;
  }
  .params-table th {
    text-align: left;
    padding: 5px 10px;
    color: var(--text-dim);
    font-weight: 400;
    border-bottom: 1px solid var(--border);
    font-size: 11px;
    text-transform: uppercase;
    letter-spacing: .4px;
  }
  .params-table td {
    padding: 6px 10px;
    border-bottom: 1px solid rgba(255,255,255,.04);
    vertical-align: top;
  }
  .param-name { color: var(--text-bright); font-weight: 500; }
  .param-ann  { color: var(--accent3); }
  .param-default { color: var(--orange); }

  /* Call graph tags */
  .callgraph-row { display: flex; flex-wrap: wrap; gap: 6px; margin-bottom: 10px; align-items: center; }
  .callgraph-label { font-size: 12px; font-weight: bold; color: var(--text-dim); min-width: 60px; }
  .call-tag {
    font-size: 11px;
    padding: 2px 8px;
    border-radius: 20px;
    cursor: pointer;
    text-decoration: none;
    transition: opacity .15s;
  }
  .call-tag:hover { opacity: .8; }
  .callee-tag { background: rgba(91,138,240,.15); color: var(--accent2); border: 1px solid rgba(91,138,240,.3); }
  .caller-tag { background: rgba(139,92,246,.15); color: var(--accent3); border: 1px solid rgba(139,92,246,.3); }
  .none-tag   { color: var(--text-dim); font-size: 11px; font-style: italic; }

  /* Class methods sub-section */
  .methods-section { margin-top: 16px; }
  .methods-label {
    font-size: 11px;
    text-transform: uppercase;
    letter-spacing: .5px;
    color: var(--text-dim);
    margin-bottom: 8px;
  }

  /* Global value */
  .global-value {
    background: var(--surface2);
    border-radius: 6px;
    padding: 6px 10px;
    font-size: 12px;
    color: var(--orange);
    margin-bottom: 10px;
    word-break: break-all;
  }
  .global-ann { font-size: 12px; color: var(--accent3); margin-bottom: 6px; }

  /* Stats bar */
  .stats-bar {
    display: flex;
    gap: 24px;
    margin-bottom: 40px;
    flex-wrap: wrap;
  }
  .stat {
    display: flex;
    flex-direction: column;
    align-items: flex-start;
  }
  .stat-num {
    font-family: 'Fraunces', serif;
    font-size: 28px;
    color: var(--accent);
    line-height: 1;
  }
  .stat-label { font-size: 11px; color: var(--text-dim); margin-top: 2px; }

  /* Lineno */
  .lineno { font-size: 10px; color: var(--text-dim); margin-left: auto; margin-top: 3px; }

  /* Scrollbar */
  ::-webkit-scrollbar { width: 6px; }
  ::-webkit-scrollbar-track { background: transparent; }
  ::-webkit-scrollbar-thumb { background: var(--border); border-radius: 3px; }

  @media (max-width: 720px) {
    #sidebar { display: none; }
    #main { padding: 24px 20px; }
  }
</style>
</head>
<body>

<nav id="sidebar">
  <div id="sidebar-header">
    <h1>{{TITLE}}</h1>
    <p>API Documentation</p>
  </div>
  <input id="search-box" type="search" placeholder="Search symbols…" autocomplete="off">
  <div id="nav-tree"></div>
</nav>

<main id="main">
  <h1 class="page-title">{{TITLE}}</h1>
  <p class="page-subtitle">Generated documentation — {{FILE_COUNT}} module(s)</p>

  <div class="stats-bar" id="stats-bar"></div>

  <div id="content"></div>
</main>

<script>
const MODULES = {{MODULES_JSON}};

// ── Stats ──────────────────────────────────────────────────────────────────
(function renderStats() {
  let totalFns = 0, totalClasses = 0, totalGlobals = 0, totalMethods = 0;
  for (const m of MODULES) {
    totalFns += m.functions.length;
    totalClasses += m.classes.length;
    totalGlobals += m.globals.length;
    for (const c of m.classes) totalMethods += c.methods.length;
  }
  const bar = document.getElementById('stats-bar');
  const items = [
    [MODULES.length, 'Modules'],
    [totalClasses, 'Classes'],
    [totalFns, 'Functions'],
    [totalMethods, 'Methods'],
    [totalGlobals, 'Globals'],
  ];
  bar.innerHTML = items.map(([n, l]) =>
    `<div class="stat"><span class="stat-num">${n}</span><span class="stat-label">${l}</span></div>`
  ).join('');
})();

// ── Helpers ────────────────────────────────────────────────────────────────
function esc(s) {
  return String(s ?? '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
}

function kindBadge(kind) {
  return `<span class="entry-kind kind-${kind}">${kind}</span>`;
}

function paramRows(params) {
  if (!params.length) return '';
  const rows = params.map(p => `
    <tr>
      <td class="param-name">${esc(p.name)}</td>
      <td class="param-ann">${p.annotation ? esc(p.annotation) : '<span style="opacity:.3">—</span>'}</td>
      <td class="param-default">${p.default ? esc(p.default) : ''}</td>
    </tr>`).join('');
  return `
    <table class="params-table">
      <thead><tr><th>Parameter</th><th>Type</th><th>Default</th></tr></thead>
      <tbody>${rows}</tbody>
    </table>`;
}

function callTags(calls, callers) {
  const calleeHtml = calls.length
    ? calls.map(c => `<a class="call-tag callee-tag" href="#${anchorId(c)}" title="Go to ${esc(c)}">${esc(c)}</a>`).join('')
    : '<span class="none-tag">none</span>';
  const callerHtml = callers.length
    ? callers.map(c => `<a class="call-tag caller-tag" href="#${anchorId(c)}" title="Go to ${esc(c)}">${esc(shortName(c))}</a>`).join('')
    : '<span class="none-tag">none</span>';
  return `
    <div class="callgraph-row">
      <span class="callgraph-label">Calls →</span>${calleeHtml}
    </div>
    <div class="callgraph-row">
      <span class="callgraph-label">← Called by</span>${callerHtml}
    </div>`;
}

function anchorId(qualName) {
  return qualName.replace(/\./g, '-');
}
function shortName(qualName) {
  const parts = qualName.split('.');
  return parts.slice(-2).join('.');
}

// ── Render function card ───────────────────────────────────────────────────
function renderFunc(f, insideClass) {
  const id = anchorId(f.qualified_name);
  const paramsStr = f.params.map(p => {
    let s = p.name;
    if (p.annotation) s += ': ' + p.annotation;
    if (p.default) s += ' = ' + p.default;
    return s;
  }).join(', ');
  const sigStr = `${esc(f.name)}(${esc(paramsStr)})`;
  const retStr = f.returns ? `→ <span>${esc(f.returns)}</span>` : '';

  return `
  <div class="entry-card" id="${id}">
    <div class="entry-header" onclick="toggleBody(this)">
      ${kindBadge(f.kind)}
      <div class="entry-sig">
        <div class="entry-name">${sigStr}</div>
        ${retStr ? `<div class="entry-return">${retStr}</div>` : ''}
      </div>
      <div class="lineno">L${f.lineno}</div>
    </div>
    <div class="entry-body">
      ${f.docstring
        ? `<div class="docstring-text">${esc(f.docstring)}</div>`
        : '<div class="no-doc">No docstring provided.</div>'}
      ${paramRows(f.params)}
      ${callTags(f.calls, f.callers)}
    </div>
  </div>`;
}

// ── Render class card ──────────────────────────────────────────────────────
function renderClass(c) {
  const id = anchorId(c.qualified_name);
  const basesStr = c.bases.length ? ` (${c.bases.map(esc).join(', ')})` : '';
  const methodsHtml = c.methods.map(m => renderFunc(m, true)).join('');

  return `
  <div class="entry-card" id="${id}">
    <div class="entry-header" onclick="toggleBody(this)">
      ${kindBadge('class')}
      <div class="entry-sig">
        <div class="entry-name">${esc(c.name)}${basesStr}</div>
      </div>
      <div class="lineno">L${c.lineno}</div>
    </div>
    <div class="entry-body">
      ${c.docstring
        ? `<div class="docstring-text">${esc(c.docstring)}</div>`
        : '<div class="no-doc">No docstring provided.</div>'}
      ${c.methods.length ? `
        <div class="methods-section">
          <div class="methods-label">Methods (${c.methods.length})</div>
          ${methodsHtml}
        </div>` : ''}
    </div>
  </div>`;
}

// ── Render global card ─────────────────────────────────────────────────────
function renderGlobal(g) {
  const id = anchorId(g.qualified_name);
  return `
  <div class="entry-card" id="${id}">
    <div class="entry-header" onclick="toggleBody(this)">
      ${kindBadge('global')}
      <div class="entry-sig">
        <div class="entry-name">${esc(g.name)}</div>
        ${g.annotation ? `<div class="entry-return"><span>${esc(g.annotation)}</span></div>` : ''}
      </div>
      <div class="lineno">L${g.lineno}</div>
    </div>
    <div class="entry-body">
      ${g.value ? `<div class="global-value">${esc(g.value)}</div>` : ''}
      ${g.docstring
        ? `<div class="docstring-text">${esc(g.docstring)}</div>`
        : '<div class="no-doc">No docstring provided.</div>'}
    </div>
  </div>`;
}

// ── Render all modules ─────────────────────────────────────────────────────
function renderAll() {
  const content = document.getElementById('content');
  content.innerHTML = MODULES.map(m => {
    const funcHtml    = m.functions.map(f => renderFunc(f, false)).join('');
    const classHtml   = m.classes.map(c => renderClass(c)).join('');
    const globalHtml  = m.globals.map(g => renderGlobal(g)).join('');

    return `
    <section class="module-section" id="mod-${anchorId(m.name)}">
      <div class="module-name">
        ${esc(m.name)}
        <span class="module-path">${esc(m.path)}</span>
      </div>
      ${m.docstring ? `<div class="module-docstring">${esc(m.docstring)}</div>` : ''}
      ${globalHtml}
      ${classHtml}
      ${funcHtml}
    </section>`;
  }).join('');
}

// ── Sidebar nav ────────────────────────────────────────────────────────────
function renderNav() {
  const tree = document.getElementById('nav-tree');
  tree.innerHTML = MODULES.map(m => {
    const items = [
      ...m.globals.map(g => ({ name: g.name, kind: 'global', id: anchorId(g.qualified_name) })),
      ...m.classes.map(c => ({ name: c.name, kind: 'class', id: anchorId(c.qualified_name) })),
      ...m.functions.map(f => ({ name: f.name, kind: 'function', id: anchorId(f.qualified_name) })),
    ];
    const itemsHtml = items.map(i => `
      <a class="nav-item" data-search="${i.name.toLowerCase()}" href="#${i.id}" onclick="navClick(event, '${i.id}')">
        <span class="kind-badge">${i.kind[0].toUpperCase()}</span>
        ${esc(i.name)}
      </a>`).join('');

    return `
    <div class="nav-module">
      <div class="nav-module-header" onclick="toggleNav(this)">
        <span class="chevron">▾</span>
        ${esc(m.name)}
      </div>
      <div class="nav-items">${itemsHtml}</div>
    </div>`;
  }).join('');
}

// ── Interactions ───────────────────────────────────────────────────────────
function toggleBody(header) {
  const body = header.nextElementSibling;
  body.classList.toggle('open');
}

function toggleNav(header) {
  header.classList.toggle('collapsed');
  const items = header.nextElementSibling;
  items.style.display = header.classList.contains('collapsed') ? 'none' : '';
}

function expandAndHighlight(id) {
  const el = document.getElementById(id);
  if (!el) return;
  // Auto-expand the card body
  const body = el.querySelector('.entry-body');
  if (body) body.classList.add('open');
  // Highlight matching nav item
  document.querySelectorAll('.nav-item').forEach(n => {
    n.classList.toggle('active', n.getAttribute('href') === '#' + id);
  });
}

function navClick(event, id) {
  event.preventDefault();
  const el = document.getElementById(id);
  if (!el) return;
  // Update URL hash without jumping
  history.pushState(null, '', '#' + id);
  el.scrollIntoView({ behavior: 'smooth', block: 'start' });
  expandAndHighlight(id);
}

// ── Search ─────────────────────────────────────────────────────────────────
document.getElementById('search-box').addEventListener('input', function() {
  const q = this.value.toLowerCase().trim();
  document.querySelectorAll('.nav-item').forEach(item => {
    const match = !q || (item.dataset.search || '').includes(q);
    item.style.display = match ? '' : 'none';
  });
});

// ── Init ───────────────────────────────────────────────────────────────────
renderAll();
renderNav();
// Honour any #hash in the URL on page load
if (location.hash) {
  const id = location.hash.slice(1);
  setTimeout(() => {
    const el = document.getElementById(id);
    if (el) {
      el.scrollIntoView({ block: 'start' });
      expandAndHighlight(id);
    }
  }, 0);
}
</script>
</body>
</html>
"""


def render_html(modules: list[ModuleDoc], title: str) -> str:
    modules_json = to_json(modules)
    html = HTML_TEMPLATE
    html = html.replace("{{TITLE}}", title)
    html = html.replace("{{FILE_COUNT}}", str(len(modules)))
    html = html.replace("{{MODULES_JSON}}", modules_json)
    return html


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Generate HTML documentation with call graphs for a Python project.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              python pydocgen.py mymodule.py
              python pydocgen.py src/ --output docs/index.html --title "MyLib Docs"
        """)
    )
    parser.add_argument("source", help="Python file or directory to document")
    parser.add_argument("-o", "--output", default="docs.html", help="Output HTML file (default: docs.html)")
    parser.add_argument("--title", default="Python Docs", help="Documentation title")
    args = parser.parse_args()

    source = Path(args.source).resolve()
    if not source.exists():
        print(f"Error: {source} does not exist.", file=sys.stderr)
        sys.exit(1)

    base_dir = source if source.is_dir() else source.parent
    files = collect_python_files(source)

    if not files:
        print("No Python files found.", file=sys.stderr)
        sys.exit(1)

    print(f"Parsing {len(files)} file(s)…")
    modules = []
    for f in files:
        print(f"  → {f.relative_to(base_dir)}")
        modules.append(parse_module(f, base_dir))

    print("Building call graph…")
    resolve_call_graph(modules)

    print("Rendering HTML…")
    html = render_html(modules, args.title)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    print(f"\n✓ Documentation written to: {out.resolve()}")
    print(f"  Open in your browser: file://{out.resolve()}")


if __name__ == "__main__":
    main()
