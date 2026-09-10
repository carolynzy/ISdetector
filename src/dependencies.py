"""
Module: src/dependencies.py
Project: ISdetector
Version: 1.1.0
Purpose: Validate external executable versions and installation-specific paths.
"""
import re
import os
import subprocess
import json
import sys

MIN_SAMTOOLS_VERSION = (1, 3)


def configured_executable(name):
    """Honor installation-specific paths; direct pip installs retain PATH discovery."""
    config_path = os.path.join(sys.prefix, 'etc', 'isdetector', 'tools.json')
    if not os.path.exists(config_path):
        return select_executable(name)
    try:
        with open(config_path) as stream:
            path = json.load(stream)['tools'][name]
        if not isinstance(path, str) or not os.path.isabs(path):
            raise ValueError('tool path must be absolute')
        if not os.path.isfile(path) or not os.access(path, os.X_OK):
            raise ValueError('saved executable is missing or not executable: {}'.format(path))
        if name == 'samtools':
            check_samtools_version(path)
        return path
    except (OSError, ValueError, KeyError, TypeError, RuntimeError) as exc:
        raise RuntimeError('Invalid tool configuration {}: {}. Rerun install.sh with a new --prefix.'.format(config_path, exc)) from exc


def executable_candidates(name):
    """Yield executable files in PATH order, ignoring duplicate paths."""
    seen = set()
    for directory in os.environ.get('PATH', os.defpath).split(os.pathsep):
        path = os.path.abspath(os.path.join(directory or os.curdir, name))
        if path not in seen and os.path.isfile(path) and os.access(path, os.X_OK):
            seen.add(path)
            yield path


def select_executable(name):
    """Select the first candidate satisfying the tool's declared requirements."""
    errors = []
    for path in executable_candidates(name):
        if name == 'samtools':
            try:
                check_samtools_version(path)
            except RuntimeError as exc:
                errors.append(str(exc))
                continue
        return path
    raise RuntimeError('No compatible {} found in PATH.\n{}'.format(name, '\n'.join(errors)))


def check_samtools_version(executable):
    minimum = '.'.join(map(str, MIN_SAMTOOLS_VERSION))
    guidance = ('ISdetector requires Samtools >= {}. Install a compatible version '
                'and ensure it is selected first in PATH.').format(minimum)
    try:
        result = subprocess.run(
            [executable, '--version'], stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, universal_newlines=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError('{} Could not query {}: {}'.format(guidance, executable, exc)) from exc

    match = re.search(r'^samtools\s+(\d+)\.(\d+)([^\s]*)', result.stdout, re.MULTILINE)
    if result.returncode != 0 or match is None:
        details = (result.stderr or result.stdout).strip()
        raise RuntimeError(
            '{} Could not determine the version of {} (exit {}).\n{}'.format(
                guidance, executable, result.returncode, details))
    version = (int(match.group(1)), int(match.group(2)))
    if version < MIN_SAMTOOLS_VERSION:
        raise RuntimeError('{} Found {} at {}.'.format(guidance, match.group(0), executable))
    return match.group(0)
