#!/usr/bin/bash

git add README.md
git commit -m "README"
git add swath-ms.ini swath-ms.R swath-ms.ipynb swath-ms.sh
git commit -m "Primary Code"
git add utils
git commit -m "utilities"
git push
