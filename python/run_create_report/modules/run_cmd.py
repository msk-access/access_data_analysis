import subprocess
import typer
from rich import print


def run_cmd(cmd):
    """Given a system command run it using subprocess

    Args:
        cmd (str): System command to be run as a string
    """
    print("\nrun_cmd:Command:\n")
    print(cmd)
    out = subprocess.Popen(
        (cmd),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
    )
    out.wait()
    stdout, stderr = out.communicate()
    if stderr is None and stdout is not None:
        if "Error" in stdout or "error" in stdout:
            print("run_cmd:stdout:\n")
            print(stdout)
            exit(1)
    else:
        print("run_cmd:stderr:\n")
        print(stderr)
    return out


def run_multiple_cmd(commands, parallel_process=None):
    """Given a system command run it using subprocess

    Args:
        cmd (list[str]): list of system commands to be run
    """
    if not parallel_process:
        parallel_process = 5
    for j in range(max(int(len(commands) / parallel_process), 1)):
        process = [
            subprocess.Popen(i, shell=True)
            for i in commands[
                j * parallel_process : min((j + 1) * parallel_process, len(commands))
            ]
        ]
    for p in process:
        p.wait()
    return ()
