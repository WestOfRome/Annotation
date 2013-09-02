#!/bin/bash 

# http://gitref.org/basic/

# http://stackoverflow.com/questions/1314950/git-get-all-commits-and-blobs-they-created
# git cat-file --batch-check < <list-of-blob-shas>
# git log --name-status

git rev-list --all --pretty=oneline

git commit -am $@ 

git rev-list --all --pretty=oneline

