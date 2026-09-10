"""
Module: src/processes.py
Project: ISdetector
Version: 1.1.0
Purpose: Run external tools with consistent diagnostics and process cleanup.
"""
import logging
import shlex
import subprocess
import threading
from collections import deque
from contextlib import contextmanager
from src.dependencies import configured_executable


def resolve_command(command):
    if command[0] in ('bwa', 'samtools', 'blastn'):
        return [configured_executable(command[0])] + list(command[1:])
    return command


def command_text(command):
    return ' '.join(shlex.quote(str(arg)) for arg in command)


def check_process(process, stderr):
    if process.returncode:
        message = "Command failed (exit {}): {}".format(
            process.returncode, command_text(process.args))
        if stderr.strip():
            message += "\n" + stderr.strip()
        logging.error(message)
        raise subprocess.CalledProcessError(
            process.returncode, process.args, stderr=stderr)


def run_checked(command):
    command = resolve_command(command)
    try:
        result = subprocess.run(command, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, universal_newlines=True)
    except OSError as exc:
        logging.error("Could not start %s: %s", command_text(command), exc)
        raise
    check_process(result, result.stderr)
    return result.stdout


@contextmanager
def streaming_process(command, stdin=None):
    """Drain stderr while streaming stdout, then check status even on read errors."""
    command = resolve_command(command)
    try:
        process = subprocess.Popen(command, stdin=stdin, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    except OSError as exc:
        logging.error("Could not start %s: %s", command_text(command), exc)
        raise
    errors = deque(maxlen=200)

    def drain():
        for line in iter(process.stderr.readline, b''):
            errors.append(line.decode('utf-8', errors='replace'))

    thread = threading.Thread(target=drain)
    thread.start()
    try:
        yield process
    except BaseException:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
        raise
    finally:
        # Closing the reader also releases producers when a downstream stage fails.
        process.stdout.close()
        process.wait()
        thread.join()
        process.stderr.close()
        check_process(process, ''.join(errors))
