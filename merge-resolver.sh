#!/bin/bash

echo "Please enter the commit hash you want to cherry-pick:"
read COMMIT_HASH

# Validate commit hash
if ! [[ $COMMIT_HASH =~ ^[0-9a-f]{40}$ ]] && ! [[ $COMMIT_HASH =~ ^[0-9a-f]{7,8}$ ]]; then
    echo "Error: Invalid commit hash format"
    exit 1
fi

# Abort existing cherry-pick if any
git cherry-pick --abort 2>/dev/null

# Start cherry-pick
echo "Starting cherry-pick for commit: $COMMIT_HASH"
git cherry-pick -n $COMMIT_HASH

# Check if cherry-pick resulted in conflicts
if git diff --name-only --diff-filter=U | grep -q .; then
    echo "Resolving conflicts..."

    # Handle each conflicted file
    git diff --name-only --diff-filter=U | while read file; do
        if [ -f "$file" ]; then
            echo "Processing: $file"

            # Get the base version
            git show :1:"$file" > "$file.base" 2>/dev/null
            # Get their version
            git show "$COMMIT_HASH":"$file" > "$file.theirs" 2>/dev/null
            # Get our version
            git show HEAD:"$file" > "$file.ours" 2>/dev/null

            # Use merge-file with correct options
            git merge-file -L "ours" -L "base" -L "theirs" --theirs --zdiff3 -p "$file.ours" "$file.base" "$file.theirs" > "$file.tmp"

            # Remove any remaining merge markers
            sed -i '' -e '/^<<<<<<< /d' \
                     -e '/^=======$/d' \
                     -e '/^>>>>>>> /d' \
                     -e '/^||||||| /d' "$file.tmp"

            # Move cleaned file to original
            mv "$file.tmp" "$file"

            # Cleanup temp files
            rm -f "$file.base" "$file.ours" "$file.theirs"

            # Stage the resolved file
            git add "$file"
        else
            git rm "$file" 2>/dev/null
        fi
    done

    # Verify all conflicts are resolved
    if git diff --check; then
        git commit -C $COMMIT_HASH
        echo "Cherry-pick completed successfully!"
    else
        echo "Error: Unresolved conflicts remain"
        exit 1
    fi
else
    git commit -C $COMMIT_HASH
    echo "Cherry-pick completed successfully!"
fi

exit 0