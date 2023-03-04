import typer
import subprocess
def run_cmd(cmd):
    """Given a system command run it using subprocess 

    Args:
        cmd (str): System command to be run as a string
    """
    typer.secho(
        f"run_cmd: run: the command line is {cmd}",
        fg=typer.colors.BRIGHT_GREEN,
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
        typer.secho(f"run_cmd: run: Read: {stdout}", 
                    fg=typer.colors.BRIGHT_GREEN,
                    )
    else:
        typer.secho(
            f"run_cmd: Could not run the command. {stderr}",
            err=True,
            fg=typer.colors.BRIGHT_RED,
            )
    return(out)