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


def get_git_tag(path="./"):
    """Get the current git tag for the repository."""
    completed_process = subprocess.run(
        ["git", "-C", os.path.abspath(path), "describe", "--tags", "--always"],
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


def get_git_reponame(path="./"):
    """Get the repository name from the remote origin URL."""
    completed_process = subprocess.run(
        ["git", "-C", os.path.abspath(path), "config", "--get", "remote.origin.url"],
        capture_output=True,
        text=True,
        cwd=os.path.dirname(os.path.abspath(__file__)),
    )
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    repo_name = completed_process.stdout.strip()
    # Remove trailing and leading information
    repo_name = os.path.splitext(os.path.basename(repo_name))[0]
    return repo_name


def submodule_info(path, outfile=None):
    """
    Get the submodule name and tag at the given path.
    Optionaly write to a file.
    """
    name = get_git_reponame(path)
    tag = get_git_tag(path)
    if outfile:
        with open(outfile, "w") as out:
            out.write(f"{name}-{tag}\n")
    return f"{name}-{tag}"


def main_info(outfile, **kwargs):
    """
    Writes the main git information to a file.
    Optinally takes paths to main dependencies and writes their basename
    """
    name = get_git_reponame()
    tag = get_git_tag()
    with open(outfile, "w") as out:
        out.write(f"{name}-{tag}\n")
        for title, path in kwargs.items():
            out.write(f"{title}: {os.path.basename(path)}\n")
