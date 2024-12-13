import subprocess
import os

def get_git_commit():
    """Get the latest commit hash"""
    completed_process = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'],
                                       capture_output=True, text=True,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    git_commit = completed_process.stdout.strip()
    return git_commit

def get_git_tag():
    """Get the current git tag for the repository."""
    completed_process = subprocess.run(['git', 'describe', '--tags'],
                                       capture_output=True, text=True,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    git_tag = completed_process.stdout.strip()
    return git_tag

def is_branch_clean():
    """Return whether current branch/commit is dirty"""
    completed_process = subprocess.run(['git', 'describe', '--tags', '--always', '--dirty'],
                                       capture_output=True, text=True,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    git_tag_or_hash = completed_process.stdout.strip()
    return 'dirty' not in git_tag_or_hash

def get_git_reponame():
    """Get the repository name from the remote origin URL."""
    completed_process = subprocess.run(['git', 'config', '--get', 'remote.origin.url'],
                                       capture_output=True, text=True,
                                       cwd=os.path.dirname(os.path.abspath(__file__)))
    if completed_process.returncode != 0:
        raise Exception(f"Error: {completed_process.stderr}")
    repo_name = completed_process.stdout.strip()
    # Remove trailing and leading information
    repo_name = os.path.basename(repo_name)[:-4]
    return repo_name