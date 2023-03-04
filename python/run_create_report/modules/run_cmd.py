import typer
import subprocess
def run_cmd(cmd):
    """Given a system command run it using subprocess 

    Args:
        cmd (str): System command to be run as a string
    """
    typer.secho(
        "run_cmd: run: the command line is %s",
        cmd.encode("unicode_escape").decode("utf-8"),
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
        typer.secho("run_cmd: run: Read: %s", 
                    stdout.decode("utf-8"), 
                    fg=typer.colors.BRIGHT_GREEN,
                    )
    else:
        typer.secho(
            "run_cmd: Could not run the command. %s",
            stderr.decode("utf-8"),
            err=True,
            fg=typer.colors.BRIGHT_RED,
            )
    return(out)