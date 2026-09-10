"""
Module: scripts/configure_installation.py
Project: ISdetector
Version: 1.1.0
Purpose: Validate a new environment and publish its private configuration and launcher.
"""
import json
import os
from pathlib import Path
import shlex
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.dependencies import check_samtools_version


def atomic_write(path, text, mode):
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(dir=str(path.parent))
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(text)
        os.chmod(temporary, mode)
        os.replace(temporary, str(path))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def configure(prefix, manager, launcher):
    prefix = Path(prefix).resolve()
    tools = {}
    for name in ('bwa', 'blastn', 'samtools'):
        path = prefix / 'bin' / name
        if not path.is_file() or not os.access(str(path), os.X_OK):
            raise RuntimeError('Missing executable: {}'.format(path))
        tools[name] = str(path)
    version = check_samtools_version(tools['samtools'])
    entry = prefix / 'bin' / 'isdetector'
    if not entry.is_file() or not os.access(str(entry), os.X_OK):
        raise RuntimeError('Missing installed command: {}'.format(entry))
    config = prefix / 'etc' / 'isdetector' / 'tools.json'
    atomic_write(config, json.dumps({'tools': tools, 'samtools_version': version}, indent=2) + '\n', 0o600)

    # Conda needs explicit streaming; Mamba run attaches streams by default.
    command = [manager, 'run']
    if Path(manager).name == 'conda':
        command.append('--no-capture-output')
    command += ['--prefix', str(prefix), str(entry)]
    text = ('#!/usr/bin/env bash\n# ISdetector managed launcher\n'
            'set -euo pipefail\n'
            'if (($# == 0)); then set -- --help; fi\n'
            'exec ' + ' '.join(shlex.quote(arg) for arg in command) + ' "$@"\n')
    atomic_write(Path(launcher), text, 0o755)
    print('Saved private tool configuration: {}'.format(config))


if __name__ == '__main__':
    configure(*sys.argv[1:])
