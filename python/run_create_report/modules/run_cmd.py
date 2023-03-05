import typer
import subprocess
import shlex

def run_cmd(cmd):
    """Given a system command run it using subprocess

    Args:
        cmd (str): System command to be run as a string
    """
    print("\nCommand:",cmd,"\n")
    typer.secho(
        f"run_cmd: command: {shlex.split(cmd)} ",
        fg=typer.colors.BRIGHT_MAGENTA,
    )
    out = subprocess.Popen(
        (cmd),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True,
    )
    out.wait()
    stdout, stderr = out.communicate()
    if stderr is None:
        typer.secho(
            f"run_cmd: {stdout} ",
            fg=typer.colors.BRIGHT_GREEN,
        )
    else:
        typer.secho(
            f"run_cmd: Could not run the command. {stderr} ",
            err=True,
            fg=typer.colors.BRIGHT_RED,
        )
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
