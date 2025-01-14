#!/bin/bash

echo "Please enter the commit hash you want to cherry-pick:"
read COMMIT_HASH

# Validate commit hash
if ! [[ $COMMIT_HASH =~ ^[0-9a-f]{40}$ ]] && ! [[ $COMMIT_HASH =~ ^[0-9a-f]{7,8}$ ]]; then
    echo "Error: Invalid commit hash format"
    exit 1
fi

# Abort existing cherry-pick
git cherry-pick --abort

# Start cherry-pick
echo "Starting cherry-pick for commit: $COMMIT_HASH"
git cherry-pick -n $COMMIT_HASH

# Accept all incoming changes including deletions
git diff --name-only --diff-filter=U | while read file; do
    if [ -f "$file" ]; then
        git checkout --theirs -- "$file"
    else
        git rm "$file"
    fi
done

# Stage and commit
git add -A
git commit -C $COMMIT_HASH

echo "Cherry-pick completed successfully!"
exit 0