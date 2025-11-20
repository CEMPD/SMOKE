#!/nas/longleaf/home/tranhuy/software/pkg/miniconda3/envs/AERMOD/bin/python

"""
smkruns.py

This script processes tcsh/csh run scripts for SMOKE emission modeling workflows.
It reads a user-supplied JSON configuration to group, annotate, and optionally validate
environment variable assignments (setenv/set) in the input script, then writes a new
script with improved structure and optional path checking.

Features:
- Groups setenv/set variable assignments by user-defined sections (prefix matching, case-insensitive)
- Preserves the original order of non-matched lines
- Annotates grouped variables under section headers (e.g., "# gsref")
- Unmatched setenv/set lines are placed under "# Other settings" in each region
- Optionally checks file paths and flags missing or duplicate assignments
- Ensures all comment lines start at column 1 (no leading whitespace before #)

Usage:
    python smkruns.py [--config PATH_TO_CONFIG] [--check-path] [--dry-run]

Arguments:
    --config/-c     Path to the JSON config file (default: smkruns_config.json in script directory)
    --check-path    Validate file paths referenced in variable assignments
    --dry-run       Show planned grouping summary without writing output

Requirements:
- Python 3.7 or newer
- No external dependencies (uses only Python standard library)

Configuration:
The config JSON must contain:
    - runscript_in: path to the input tcsh/csh script
    - runscript_out: path to the output script
    - sections: dictionary mapping section names to lists of variable prefixes
    - directory_definitions: path to a definitions script for environment variables

See smkruns_config.json for a template configuration.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
from typing import Dict, List, Tuple
import pathlib
import shlex


def load_config(path: str) -> Dict:
    with open(path, 'r', encoding='utf-8') as fh:
        return json.load(fh)


def extract_var_from_line(line: str) -> str | None:
    """Return the variable name assigned by `setenv` or `set` in a tcsh line.

    Examples matched:
      setenv GSREFTMP_A /some/path
      set VAR = ( ... )
    Only the simple `setenv NAME` and `set NAME` forms are detected; returns None
    if no variable name is found.
    """
    m = re.match(r"^\s*(?:setenv|set)\s+([A-Za-z0-9_]+)\b", line)
    if m:
        return m.group(1)
    return None


def parse_tcsh_env_vars(path: str, vars_of_interest: List[str]) -> Dict[str, str]:
    """Parse a tcsh/csh script for setenv or set assignments of selected variables.

    Returns a mapping of variable names to their (unquoted) values. If a variable
    is assigned multiple times, the last assignment wins.
    """
    out: Dict[str, str] = {}
    if not os.path.exists(path):
        return out
    pattern_setenv = re.compile(r"^\s*setenv\s+([A-Za-z0-9_]+)\s+(.+?)\s*$")
    pattern_set = re.compile(r"^\s*set\s+([A-Za-z0-9_]+)\s*=\s*(.+?)\s*$")
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            m1 = pattern_setenv.match(line)
            m2 = pattern_set.match(line) if not m1 else None
            if not m1 and not m2:
                continue
            var = (m1 or m2).group(1)
            if var not in vars_of_interest:
                continue
            val_part = (m1 or m2).group(2).strip()
            # Remove inline comments (## or #) not inside quotes using shlex split
            try:
                tokens = shlex.split(val_part, posix=True)
                if tokens:
                    value = tokens[0]
                else:
                    value = ''
            except Exception:
                # Fallback: strip surrounding quotes manually
                value = val_part.strip('"')
            out[var] = value
    return out

def parse_tcsh_all_env_vars(path: str) -> Dict[str, str]:
    """Parse all setenv/set assignments in a tcsh script, returning last value per var."""
    out: Dict[str, str] = {}
    if not os.path.exists(path):
        return out
    pattern_setenv = re.compile(r"^\s*setenv\s+([A-Za-z0-9_]+)\s+(.+?)\s*$")
    pattern_set = re.compile(r"^\s*set\s+([A-Za-z0-9_]+)\s*=\s*(.+?)\s*$")
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            m1 = pattern_setenv.match(line)
            m2 = pattern_set.match(line) if not m1 else None
            if not m1 and not m2:
                continue
            var = (m1 or m2).group(1)
            val_part = (m1 or m2).group(2).strip()
            try:
                tokens = shlex.split(val_part, posix=True)
                value = tokens[0] if tokens else ''
            except Exception:
                value = val_part.strip('"')
            out[var] = value
    return out


def group_region(lines: List[str], sections: Dict[str, List[str]], start: int, end: int) -> Tuple[Dict[str, List[Tuple[int, str]]], List[Tuple[int, str]], set[int]]:
    """Group lines within the [start, end) slice (indices relative to full `lines`).

    Returns grouped sections, an "other" collection for unmatched lines, and the
    set of line indices to remove from the original content.
    """
    grouped: Dict[str, List[Tuple[int, str]]] = {k: [] for k in sections.keys()}
    other: List[Tuple[int, str]] = []
    remove_idxs: set[int] = set()

    if start < 0 or end <= start:
        return grouped, other, remove_idxs

    for idx in range(start, min(end, len(lines))):
        line = lines[idx]
        var = extract_var_from_line(line)
        if var is None:
            continue
        matched = False
        for sec, patterns in sections.items():
            for p in patterns:
                if var == p or var.startswith(p):
                    print(f'[DEBUG] Grouping variable {var} (line {idx}) under section {sec}')
                    grouped[sec].append((idx, line))
                    remove_idxs.add(idx)
                    matched = True
                    break
            if matched:
                break
        if not matched:
            print(f'[DEBUG] Variable {var} (line {idx}) not grouped, goes to Other')
            other.append((idx, line))
            remove_idxs.add(idx)
    return grouped, other, remove_idxs


def render_grouped(grouped: Dict[str, List[Tuple[int, str]]], other: List[Tuple[int, str]]) -> List[str]:
    rendered: List[str] = []
    for sec, arr in grouped.items():
        if arr:
            # Interpret '\n' in section names as a real line break, prefix each with '#'
            for sline in sec.split('\\n'):
                sline = sline.rstrip('\r\n')
                if sline.strip():
                    rendered.append(f"# {sline.strip()}\n")
            for _idx, line in sorted(arr, key=lambda x: x[0]):
                rendered.append(line)
    if other:
        rendered.append("# Other settings\n")
        for _idx, line in sorted(other, key=lambda x: x[0]):
            rendered.append(line)
    return rendered


def check_region_abbrev_in_griddesc(defs_path: str, lines: List[str]) -> Tuple[str | None, str | None, bool | None]:
    """Check whether REGION_ABBREV from directory_definitions exists in the GRIDDESC file.

    Returns (region_abbrev, griddesc_path, found) where found is True/False or None if unknown.
    """
    try:
        defs_vars = parse_tcsh_all_env_vars(defs_path)
    except Exception:
        defs_vars = {}
    region_abbrev = defs_vars.get('REGION_ABBREV')

    # Find GRIDDESC assignment in the runscript (first occurrence)
    griddesc_path = None
    griddesc_pattern = re.compile(r"^\s*setenv\s+GRIDDESC\s+(.+?)\s*$")
    for line in lines:
        m = griddesc_pattern.match(line)
        if m:
            griddesc_path = m.group(1).strip().strip('"')
            break

    # Expand variables in griddesc_path using defs_vars
    if griddesc_path and defs_vars:
        var_pattern = re.compile(r"\$\{?([A-Za-z0-9_]+)\}?")
        def expand_vars(val, env):
            result = val
            for _ in range(10):
                matches = var_pattern.findall(result)
                if not matches:
                    break
                for v in set(matches):
                    if v in env:
                        result = re.sub(fr"\${{?{v}}}?", env[v], result)
            return result
        try:
            griddesc_path = expand_vars(griddesc_path, defs_vars)
        except Exception:
            pass

    # If we have both, check
    if region_abbrev and griddesc_path and os.path.exists(griddesc_path):
        found = False
        try:
            with open(griddesc_path, 'r', encoding='utf-8', errors='ignore') as gfh:
                for gl in gfh:
                    if region_abbrev in gl:
                        found = True
                        break
        except Exception:
            return region_abbrev, griddesc_path, None
        return region_abbrev, griddesc_path, found

    # Cases where we couldn't determine
    return region_abbrev, griddesc_path, None


def get_griddesc_annotation(defs_path: str, lines: list[str], check_path: bool, dry_run: bool) -> str | None:
    """Return annotation string for GRIDDESC setenv line based on REGION_ABBREV validity."""
    if not check_path or dry_run or defs_path is None or lines is None:
        return None
    region_abbrev, griddesc_path, found = check_region_abbrev_in_griddesc(defs_path, lines)
    if region_abbrev and griddesc_path:
        if found is True:
            #return f"REGION_ABBREV '{region_abbrev}' is VALID in GRIDDESC"
            # Mute annotation
            return None
        elif found is False:
            return f"REGION_ABBREV '{region_abbrev}' is NOT FOUND in GRIDDESC"
        elif found is None:
            return f"Could not read GRIDDESC file: {griddesc_path}"
    elif region_abbrev:
        return f"Could not determine GRIDDESC path from script or file missing."
    return None


def main(argv: List[str] | None = None) -> int:
    # --- Load config ---
    defs_path = None
    lines = None

    p = argparse.ArgumentParser(description="Group lines in a tcsh runscript by prefixes from a JSON config.")
    p.add_argument('--config', '-c', default=os.path.join(os.path.dirname(__file__), 'smkruns_config.json'),
                   help='Path to smkruns_config.json')
    p.add_argument('--dry-run', action='store_true', help='Parse and show planned grouping summary without writing output file.')
    p.add_argument('--check-path', dest='check_path', action='store_true', help='Validate paths of setenv lines that reference $GE_DAT (e.g., GRIDDESC).')
    p.add_argument('--no-check-path', dest='check_path', action='store_false', help='Disable path checking.')
    p.set_defaults(check_path=None)
    args = p.parse_args(argv)

    # Load config and assign variables only once after parsing args
    cfg_path = os.path.abspath(args.config)
    if not os.path.exists(cfg_path):
        print(f"Config file not found: {cfg_path}", file=sys.stderr)
        return 2

    cfg = load_config(cfg_path)
    run_in = cfg.get('runscript_in')
    run_out = cfg.get('runscript_out')
    sections = cfg.get('sections')
    defs_path = cfg.get('directory_definitions')
    config_check_path = cfg.get('check_path', False)
    # If CLI flag is not set, use config value
    if args.check_path is None:
        args.check_path = config_check_path

    if not run_in or not run_out or not sections or not defs_path:
        print(f"Config JSON must contain 'runscript_in', 'runscript_out', 'sections', and 'directory_definitions'", file=sys.stderr)
        return 3

    run_in = os.path.expanduser(run_in)
    run_out = os.path.expanduser(run_out)
    defs_path = os.path.expanduser(defs_path)

    if not os.path.exists(run_in):
        print(f"Input runscript not found: {run_in}", file=sys.stderr)
        print("No changes made. If you want to run this script without the input file present, create the input or adjust the config.")
        return 4

    if not os.path.exists(defs_path):
        print(f"Config error: directory_definitions file not found: {defs_path}", file=sys.stderr)
        return 5

    with open(run_in, 'r', encoding='utf-8') as fh:
        lines = fh.readlines()

    # --- REGION_ABBREV and GRIDDESC check ---
    if defs_path is not None and lines is not None:
        region_abbrev, griddesc_path, found = check_region_abbrev_in_griddesc(defs_path, lines)
        if region_abbrev and griddesc_path and found is True:
            print(f"[REGION_ABBREV check] REGION_ABBREV '{region_abbrev}' found in GRIDDESC: {griddesc_path}")
        elif region_abbrev and griddesc_path and found is False:
            print(f"[REGION_ABBREV check][WARNING] REGION_ABBREV '{region_abbrev}' NOT found in GRIDDESC: {griddesc_path}")
        elif region_abbrev and griddesc_path and found is None:
            print(f"[REGION_ABBREV check][WARNING] Could not read GRIDDESC file: {griddesc_path}")
        elif region_abbrev:
            print(f"[REGION_ABBREV check][WARNING] Could not determine GRIDDESC path from script or file missing.")
    # --- END REGION_ABBREV and GRIDDESC check ---

    griddesc_annotation = get_griddesc_annotation(defs_path, lines, args.check_path, args.dry_run)

    env_from_defs: Dict[str,str] = parse_tcsh_all_env_vars(defs_path)

    # Optional: path validation for lines referencing GE_DAT
    if args.check_path:
        pattern = re.compile(r"^\s*setenv\s+([A-Za-z0-9_]+)\s+(.+?)\s*$")
        ok_paths: List[str] = []
        missing_paths: List[str] = []
        skipped: List[str] = []
        checked: List[Tuple[str,str]] = []
        classifications: Dict[int, str] = {}
        # --- Duplicate detection logic ---
        var_assignments: Dict[str, List[int]] = {}
        setenv_set_pattern = re.compile(r"^\s*(setenv|set)\s+([A-Za-z0-9_]+)\b", re.IGNORECASE)
        for idx, raw in enumerate(lines):
            # Ignore lines that are fully commented out (start with # or ## after optional whitespace)
            stripped = raw.lstrip()
            if stripped.startswith('#'):
                continue
            m = setenv_set_pattern.match(raw)
            if m:
                var = m.group(2).upper()  # Case-insensitive duplicate detection
                var_assignments.setdefault(var, []).append(idx)
        duplicates = {var: idxs for var, idxs in var_assignments.items() if len(idxs) > 1}
        duplicate_line_map = {}  # Map: line idx -> (var, first_idx)
        if duplicates:
            print("[check-path][WARNING] Duplicate variable assignments detected:")
            for var, idxs in duplicates.items():
                print(f"  Variable '{var}' assigned on lines: {', '.join(str(i+1) for i in idxs)}")
                # For each duplicate after the first, record for annotation
                first_idx = idxs[0]
                for dup_idx in idxs[1:]:
                    duplicate_line_map[dup_idx] = (var, first_idx)
        # --- End duplicate detection ---
        if not env_from_defs:
            print('[check-path] No environment variables parsed from directory_definitions; skipping.')
        else:
            var_pattern = re.compile(r"\$\{?([A-Za-z0-9_]+)\}?")
            def recursive_expand(val: str, env: Dict[str, str], max_depth: int = 10) -> str:
                result = val
                for _ in range(max_depth):
                    matches = var_pattern.findall(result)
                    if not matches:
                        break
                    for v in set(matches):
                        if v in env:
                            # Replace all ${VAR} and $VAR
                            result = re.sub(fr"\${{?{v}}}?", env[v], result)
                return result

            for idx, raw in enumerate(lines):
                m = pattern.match(raw)
                if not m:
                    continue
                var, val = m.group(1), m.group(2).strip()
                # Recursively expand all env variables
                expanded = recursive_expand(val, env_from_defs)
                expanded_clean = expanded.strip().strip('"').strip("'")
                if ' ' in expanded_clean and not (expanded_clean.startswith('/') and os.path.exists(expanded_clean)):
                    expanded_clean = expanded_clean.split()[0]
                if re.match(r"^(\/|\.\.?\/)", expanded_clean):
                    exists = os.path.exists(expanded_clean)
                    checked.append((var, expanded_clean))
                    if exists:
                        ok_paths.append(f"{var}={expanded_clean}")
                        classifications[idx] = 'VALID'
                    else:
                        missing_paths.append(f"{var}={expanded_clean}")
                        classifications[idx] = 'MISSING'
                else:
                    skipped.append(f"{var}={expanded_clean}")
                    classifications[idx] = 'SKIPPED'
            print('[check-path] Summary of variable-referenced setenv lines:')
            for var, path in checked:
                print(f"  {var}: {path} -> {'OK' if os.path.exists(path) else 'MISSING'}")
            print(f"[check-path] OK: {len(ok_paths)}  MISSING: {len(missing_paths)}  SKIPPED(non-path): {len(skipped)}")
            if missing_paths:
                print('[check-path] Missing paths:')
                for item in missing_paths:
                    print(f"    {item}")
            if skipped:
                print('[check-path] Skipped (not treated as file system paths):')
                for item in skipped:
                    print(f"    {item}")


    # Define region header regex patterns for robust matching
    region_patterns = [
        ("inputs for all sectors", re.compile(r"^\s*#+?\s*inputs for all sectors", re.IGNORECASE)),
        ("inputs specific to this sector", re.compile(r"^\s*#+?\s*inputs specific to this sector", re.IGNORECASE)),
        ("parameters for all sectors", re.compile(r"^\s*#+?\s*parameters for all sectors", re.IGNORECASE)),
        ("parameters specific to this sector", re.compile(r"^\s*#+?\s*parameters specific to this sector", re.IGNORECASE)),
    ]
    # Find indices of region headers using regex
    region_indices = []
    for i, line in enumerate(lines):
        for header, pattern in region_patterns:
            if pattern.search(line):
                region_indices.append((header, i))
    # Sort by index
    region_indices.sort(key=lambda x: x[1])
    # Build region boundaries: (header, start_idx, end_idx)
    region_bounds = []
    for idx, (header, start) in enumerate(region_indices):
        end = region_indices[idx+1][1] if idx+1 < len(region_indices) else len(lines)
        region_bounds.append((header, start+1, end))

    # For each region, group variables
    grouped_regions = []
    remove_idxs = set()
    for header, start, end in region_bounds:
        grouped, other, remove = group_region(lines, sections, start, end)
        grouped_regions.append((header, start-1, grouped, other))  # start-1 is the header line idx
        remove_idxs |= remove

    # Build output lines, inserting grouped variables after each header
    out_lines = []
    input_to_output_idx = {}
    output_line_num = 0
    region_map = {header: (header, header_idx, grouped, other) for (header, header_idx, grouped, other) in grouped_regions}
    region_header_idxs = {header_idx: header for (header, header_idx, grouped, other) in grouped_regions}
    idx = 0
    while idx < len(lines):
        if idx in remove_idxs:
            idx += 1
            continue
        out_lines.append(lines[idx])
        input_to_output_idx[idx] = output_line_num
        output_line_num += 1
        # If this is a region header, insert grouped lines after it
        if idx in region_header_idxs:
            header = region_header_idxs[idx]
            grouped, other = region_map[header][2], region_map[header][3]
            grouped_lines = render_grouped(grouped, other)
            out_lines.extend(grouped_lines)
            output_line_num += len(grouped_lines)
        idx += 1

    output_lines = out_lines

    # --- Replace 'source ../directory_definitions.csh' with config value in output_lines ---
    if defs_path is not None:
        new_output_lines = []
        for line in output_lines:
            line_no_comment = re.split(r'#|;', line, 1)[0].strip()
            if re.match(r"^source\s+[^#;\n]*directory_definitions\.csh(\s|$)", line_no_comment):
                new_output_lines.append(f"source {defs_path}\n")
            else:
                new_output_lines.append(line)
        output_lines = new_output_lines

    # --- Post-process output_lines for --check-path and duplicate annotation ---
    if args.check_path:
        # Build duplicate map (case-insensitive)
        var_seen = {}
        setenv_set_pattern = re.compile(r"^\s*(setenv|set)\s+([A-Za-z0-9_]+)\b", re.IGNORECASE)
        duplicate_lines = set()
        for idx, line in enumerate(output_lines):
            m = setenv_set_pattern.match(line)
            if m:
                var = m.group(2).upper()
                if var in var_seen:
                    duplicate_lines.add(idx)
                else:
                    var_seen[var] = idx
        # Path check and annotation
        env_from_defs: Dict[str,str] = parse_tcsh_all_env_vars(defs_path)
        var_pattern = re.compile(r"\$\{?([A-Za-z0-9_]+)\}?")
        def recursive_expand(val: str, env: Dict[str, str], max_depth: int = 10) -> str:
            result = val
            for _ in range(max_depth):
                matches = var_pattern.findall(result)
                if not matches:
                    break
                for v in set(matches):
                    if v in env:
                        result = re.sub(fr"\${{?{v}}}?", env[v], result)
            return result
        pattern = re.compile(r"^\s*setenv\s+([A-Za-z0-9_]+)\s+(.+?)\s*$")
        for idx, line in enumerate(output_lines):
            orig_line = line.rstrip('\n')
            # Remove any previous VALID PATH or MISSING PATH or DUPLICATED SETTING or REGION_ABBREV annotation from the comment
            if '#' in orig_line and not orig_line.lstrip().startswith('#'):
                base, comment = orig_line.split('#', 1)
                # Remove known annotation phrases
                comment = re.sub(r'(;?\s*(VALID PATH|MISSING PATH|DUPLICATED SETTING|REGION_ABBREV[^;#\n]*))*', '', comment)
                orig_line = base.rstrip()
                if comment.strip():
                    orig_line += ' # ' + comment.strip()
            m = pattern.match(orig_line)
            annotation = ''
            if m:
                var, val = m.group(1), m.group(2).strip()
                expanded = recursive_expand(val, env_from_defs)
                expanded_clean = expanded.strip().strip('"').strip("'")
                if ' ' in expanded_clean and not (expanded_clean.startswith('/') and os.path.exists(expanded_clean)):
                    expanded_clean = expanded_clean.split()[0]
                if re.match(r"^(\/|\.\.?\/)", expanded_clean):
                    exists = os.path.exists(expanded_clean)
                    if exists:
                        annotation = 'VALID PATH'
                    else:
                        annotation = 'MISSING PATH'
                else:
                    annotation = ''
            if idx in duplicate_lines:
                if annotation:
                    annotation += ' ; DUPLICATED SETTING'
                else:
                    annotation = 'DUPLICATED SETTING'
            if annotation:
                line = orig_line.rstrip()
                if '#' in line:
                    line += f' ; {annotation}'
                else:
                    line += f' # {annotation}'
                line += '\n'
                output_lines[idx] = line
            else:
                output_lines[idx] = orig_line + '\n'
    else:
        # If --check-path is not specified, remove only known path/region annotations from setenv/set lines, but preserve all other comments (e.g., error messages)
        for idx, line in enumerate(output_lines):
            orig_line = line.rstrip('\n')
            if '#' in orig_line:
                print(f'[DEBUG] Orignal line with #-->{orig_line}')
                base, comment = orig_line.split('#', 1)
                # Remove only known annotation phrases, but preserve other comments
                # Remove all instances of these phrases, but keep other comment content
                cleaned_comment = re.sub(r'(;?\s*(VALID PATH|MISSING PATH|DUPLICATED SETTING|REGION_ABBREV[^;#\n]*))*', '', comment)
                # If the comment is now empty, remove the trailing #, else preserve the rest
                orig_line = base.rstrip()
                if cleaned_comment.strip():
                    orig_line += '# ' + cleaned_comment.strip()
            output_lines[idx] = orig_line + '\n'

    # Ensure output directory exists
    out_dir = os.path.dirname(run_out)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    if args.dry_run:
        # Print a summary of counts per section for quick inspection using grouped_regions
        summary: Dict[str, int] = {}
        for header, header_idx, grouped, other in grouped_regions:
            # Normalize header as key prefix
            key_prefix = header.replace(' ', '_')
            for sec, arr in grouped.items():
                summary[f"{key_prefix}:{sec}"] = len(arr)
            summary[f"{key_prefix}:other"] = len(other)
        print("Dry-run summary (lines grouped per section):")
        for k in sorted(summary):
            print(f"  {k}: {summary[k]}")
        print("(No output file written; remove --dry-run to write.)")
    else:
        # If annotation is needed, append to GRIDDESC setenv line in output_lines, preserving any path annotation
        if griddesc_annotation:
            pattern = re.compile(r"^(\s*setenv\s+GRIDDESC\s+.+?)(\s*#.*)?$", re.IGNORECASE)
            for i, line in enumerate(output_lines):
                m = pattern.match(line)
                if m:
                    base = m.group(1)
                    comment = m.group(2) or ''
                    # If comment already contains REGION_ABBREV annotation, remove it
                    comment = re.sub(r";?\s*REGION_ABBREV.*", "", comment)
                    # Append REGION_ABBREV annotation after any existing comment
                    if comment.strip():
                        output_lines[i] = f"{base}{comment}; {griddesc_annotation.strip()}\n"
                    else:
                        output_lines[i] = f"{base} # {griddesc_annotation}\n"
                    break
        with open(run_out, 'w', encoding='utf-8') as fh:
            fh.writelines(output_lines)
        print(f"Wrote grouped runscript to: {run_out}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
