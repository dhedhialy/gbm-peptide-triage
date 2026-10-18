#!/usr/bin/env bash
# Publish dist/ to the gh-pages branch that GitHub Pages serves.
# Run it after `python3 scripts/train.py` regenerates public/model.json.
set -euo pipefail
cd "$(dirname "$0")/.."

remote=$(git remote get-url origin)
repo=$(gh repo view --json nameWithOwner -q .nameWithOwner)
name=$(git config user.name)
email=$(git config user.email)

yarn build

tmp=$(mktemp -d)
cp -r dist/. "$tmp/"
touch "$tmp/.nojekyll"   # ship the asset filenames as-is, no Jekyll pass
git -C "$tmp" init -q -b gh-pages
git -C "$tmp" add -A
git -C "$tmp" -c user.name="$name" -c user.email="$email" \
  commit -qm "Deploy $(date -u +%Y-%m-%dT%H:%M:%SZ)"
git -C "$tmp" push -f "$remote" gh-pages
rm -rf "$tmp"

echo "pushed gh-pages; serving at https://${repo%%/*}.github.io/${repo##*/}/"
