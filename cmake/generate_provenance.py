#!/usr/bin/env python3
"""Capture Git source status at build time; write header only when changed."""
import argparse
import pathlib
import subprocess

def git(repo, *args):
    try:
        result=subprocess.run(['git','-C',str(repo),*args], capture_output=True,
                              text=True,check=False)
        return result.stdout.strip() if result.returncode==0 else None
    except OSError:
        return None

def status(repo):
    revision=git(repo,'rev-parse','--verify','HEAD') or 'unknown'
    # Tracked changes only: untracked scratch artifacts do not affect compiled source.
    # A missing repository is conservatively treated as dirty/unknown.
    changed=git(repo,'status','--porcelain','--untracked-files=no')
    return revision, changed is None or bool(changed)

def literal(value):
    # Python JSON string literals are valid C++ string literals for these values.
    import json
    return json.dumps(value, ensure_ascii=True)

p=argparse.ArgumentParser()
p.add_argument('--statefinder',required=True)
p.add_argument('--solver',required=True)
p.add_argument('--output',required=True)
p.add_argument('--build-type',default='unknown')
p.add_argument('--compiler',default='unknown')
a=p.parse_args()
sr,sd=status(a.statefinder)
fr,fd=status(a.solver)
content='''// Generated at build time. Do not edit.
#pragma once
#include "tpsf/engineering_report.hpp"
namespace tpsf {
inline EngineeringProvenance compiled_engineering_provenance() {
    EngineeringProvenance p;
'''
content += '    p.statefinder_revision = ' + literal(sr) + ';\n'
content += '    p.solver_revision = ' + literal(fr) + ';\n'
content += '    p.statefinder_dirty = ' + ('true' if sd else 'false') + ';\n'
content += '    p.solver_dirty = ' + ('true' if fd else 'false') + ';\n'
content += '    p.build_type = ' + literal(a.build_type or 'unknown') + ';\n'
content += '    p.compiler = ' + literal(a.compiler or 'unknown') + ';\n'
content += '    return p;\n}\n} // namespace tpsf\n'
out=pathlib.Path(a.output)
out.parent.mkdir(parents=True,exist_ok=True)
if not out.exists() or out.read_text(encoding='utf-8') != content:
    out.write_text(content,encoding='utf-8')
