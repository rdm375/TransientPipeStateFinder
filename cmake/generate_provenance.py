#!/usr/bin/env python3
"""Capture Git source status at build time; write header only when changed."""
import argparse
import hashlib
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

def fingerprint(repo):
    """Hash source candidates by relative path and bytes, independent of mtimes.

    Includes tracked and nonignored untracked source files in known source
    directories, plus root build manifests. This is a source-input fingerprint,
    not proof of the exact set of translation units linked by CMake.
    """
    repo=pathlib.Path(repo).resolve()
    try:
        result=subprocess.run(
            ['git','-C',str(repo),'ls-files','--cached','--others',
             '--exclude-standard','-z'], capture_output=True, check=True)
    except (OSError, subprocess.CalledProcessError):
        return 'unknown'
    allowed={'.cpp','.cc','.cxx','.c','.cu','.cuh','.hpp','.hh','.h',
             '.hxx','.f','.f90','.f95','.f03','.f08','.for','.py','.cmake'}
    directories={'src','include','apps','cmake','fortran','cuda','kernels',
                 'tests','lib','modules'}
    manifests={'CMakeLists.txt','CMakePresets.json','Makefile','meson.build'}
    digest=hashlib.sha256()
    paths=sorted(set(result.stdout.split(b'\0'))-{b''})
    for raw in paths:
        rel=pathlib.Path(raw.decode('utf-8','surrogateescape'))
        if rel.is_absolute() or '..' in rel.parts:
            continue
        if not (rel.name in manifests or
                (rel.parts[0] in directories and
                 (rel.suffix.lower() in allowed or rel.name in manifests))):
            continue
        path=repo/rel
        if not path.is_file() or path.is_symlink():
            continue
        content=path.read_bytes()
        digest.update(len(raw).to_bytes(8,'big'))
        digest.update(raw)
        digest.update(len(content).to_bytes(8,'big'))
        digest.update(content)
    return digest.hexdigest()

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
sf=fingerprint(a.statefinder)
ff=fingerprint(a.solver)
content='''// Generated at build time. Do not edit.
#pragma once
#include "tpsf/engineering_report.hpp"
namespace tpsf {
inline EngineeringProvenance compiled_engineering_provenance() {
    EngineeringProvenance p;
'''
content += '    p.statefinder_revision = ' + literal(sr) + ';\n'
content += '    p.solver_revision = ' + literal(fr) + ';\n'
content += '    p.statefinder_source_sha256 = ' + literal(sf) + ';\n'
content += '    p.solver_source_sha256 = ' + literal(ff) + ';\n'
content += '    p.statefinder_dirty = ' + ('true' if sd else 'false') + ';\n'
content += '    p.solver_dirty = ' + ('true' if fd else 'false') + ';\n'
content += '    p.build_type = ' + literal(a.build_type or 'unknown') + ';\n'
content += '    p.compiler = ' + literal(a.compiler or 'unknown') + ';\n'
content += '    return p;\n}\n} // namespace tpsf\n'
out=pathlib.Path(a.output)
out.parent.mkdir(parents=True,exist_ok=True)
if not out.exists() or out.read_text(encoding='utf-8') != content:
    out.write_text(content,encoding='utf-8')
