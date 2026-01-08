import subprocess
import os


def get_git_commit():
    """Get the latest commit hash"""
    completed_process = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.abspath(__file__)),
    )
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    git_commit = completed_process.stdout.strip()
    return git_commit


def get_git_tag(path=None):
    """Get the current git tag for the repository."""
    if path is None:
        cmd = ["git", "describe", "--tags", "--always"]
    else:
        cmd = ["git", "-C", os.path.abspath(path), "describe", "--tags", "--always"]

    completed_process = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.abspath(__file__)),
    )
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    git_tag = completed_process.stdout.strip()
    return git_tag


def is_branch_clean():
    """Return whether current branch/commit is dirty"""
    completed_process = subprocess.run(
        ["git", "describe", "--tags", "--always", "--dirty"],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.abspath(__file__)),
    )
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    git_tag_or_hash = completed_process.stdout.strip()
    return "dirty" not in git_tag_or_hash


def get_git_reponame(path=None):
    """
    Get the repository name from the remote origin URL.

    If path is None, use the repo containing this file.
    If path is given, treat it as the repo root.
    """
    if path is None:
        cmd = ["git", "config", "--get", "remote.origin.url"]
    else:
        cmd = [
            "git",
            "-C",
            os.path.abspath(path),
            "config",
            "--get",
            "remote.origin.url",
        ]

    completed_process = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.abspath(__file__)),
    )

    if completed_process.returncode != 0:
        raise Exception(
            f"Error running git in get_git_reponame "
            f"(path={path!r}): {completed_process.stderr}"
        )

    repo_url = completed_process.stdout.strip()
    repo_name = os.path.splitext(os.path.basename(repo_url))[0]
    return repo_name


def submodule_info(path, outfile=None):
    """
    Get the submodule name and tag at the given path.
    Optionally write to a file.
    """
    name = get_git_reponame(path)
    tag = get_git_tag(path)
    if outfile:
        with open(outfile, "a") as out:
            out.write(f"{name}-{tag}\n")
    return f"{name}-{tag}"


def main_info(outfile, **kwargs):
    """
    Writes the main git information to a file.
    Optionally takes paths to main dependencies and writes their basename
    """
    name = get_git_reponame()
    tag = get_git_tag()
    with open(outfile, "a") as out:
        out.write(f"{name}-{tag}\n")
        for title, path in kwargs.items():
            out.write(f"{title}: {os.path.basename(path)}\n")
