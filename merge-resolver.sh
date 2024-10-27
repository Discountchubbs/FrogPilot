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

# Define excluded patterns
EXCLUDED_PATTERNS=(
    "hyundai/"
    "controlsd"
    "radard"
    "hyundairadartracks"
    "hkgtuning"
    "hattrick"
    "brakesignal"
)

# Process each conflicted file
git diff --name-only --diff-filter=U | while read file; do
    if [ -f "$file" ]; then
        echo "Processing: $file"

        # Get the three versions of the file
        git show :1:"$file" > "$file.base"    # common ancestor
        git show :2:"$file" > "$file.ours"    # our changes
        git show :3:"$file" > "$file.theirs"  # their changes

        # Check if file contains excluded patterns
        CONTAINS_EXCLUDED=0
        for pattern in "${EXCLUDED_PATTERNS[@]}"; do
            if grep -q "$pattern" "$file.ours"; then
                CONTAINS_EXCLUDED=1
                break
            fi
        done

        if [ $CONTAINS_EXCLUDED -eq 1 ]; then
            echo "File contains excluded patterns, performing selective merge: $file"

            # Create patch without excluded patterns
            for pattern in "${EXCLUDED_PATTERNS[@]}"; do
                # Store lines containing excluded patterns
                grep -n "$pattern" "$file.ours" | cut -d: -f1 > "$file.excluded_lines"

                # Create patch excluding protected lines
                git diff "$file.base" "$file.theirs" | \
                while IFS= read -r line; do
                    LINE_NUM=$(echo "$line" | grep -o '^@@ .* @@' | cut -d' ' -f3 | cut -d, -f1)
                    if ! grep -q "^$LINE_NUM$" "$file.excluded_lines" 2>/dev/null; then
                        echo "$line"
                    fi
                done > "$file.safe.patch"

                # Apply safe patch
                patch "$file.ours" "$file.safe.patch" -o "$file"
            done
        else
            echo "Accepting all changes for: $file"
            cp "$file.theirs" "$file"
        fi

        # Cleanup temporary files
        rm -f "$file.base" "$file.ours" "$file.theirs" "$file.excluded_lines" "$file.safe.patch"
    else
        echo "Removing file: $file"
        git rm "$file"
    fi
done

# Stage and commit
git add -A
git commit -C $COMMIT_HASH

echo "Cherry-pick completed successfully!"
exit 0