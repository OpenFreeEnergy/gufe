# A simple script to handle raising or closing a github issue
# Is controlled by the following preset environment variables:
#  - GITHUB_TOKEN: A github API token
#  - CI_OUTCOME: The outcome status of a Github action step
#  - TITLE: A string title which the issue will have.

import os

from github import Github

if __name__ == "__main__":
    git = Github(os.environ["GITHUB_TOKEN"])
    status = os.environ["CI_OUTCOME"]
    repo = git.get_repo("OpenFreeEnergy/gufe")
    title = os.environ["TITLE"]

    target_issue = None
    for issue in repo.get_issues():
        if issue.title == title:
            target_issue = issue

    # Close any issues with given title if CI returned green
    if status == "success":
        if target_issue is not None:
            target_issue.edit(state="closed")
    else:
        # Otherwise raise an issue
        if target_issue is None:
            repo.create_issue(title=title, body="")
