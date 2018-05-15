#!/bin/sh

# Note: relative path is better, so that we can move the repository if we want.
hookdir="`dirname "$0"`"

# Note: We make the assumption that this script is located in .githooks/ in the
# root folder of the repository.
gitdir="$hookdir/../.git"
if test -f "$gitdir"; then
	gitdir="`grep -e '^gitdir: ' "$gitdir"`"
	gitdir="${gitdir##gitdir: }"
	if test "${gitdir:0:1}" != "/"; then # not an absolute path
		gitdir="../$gitdir"
	fi
fi

# Try to canonicalize. Note that Mac OS X does not support "readlink -f"
gitdir="`readlink -f "$gitdir" || echo "$gitdir"`"

echo "Entering $gitdir/hooks..."
cd "$gitdir/hooks"
echo "Entered `pwd`..."

for i in ../../.githooks/*; do
	# ignore myself
	if test "`basename "$i"`" != "`basename "$0"`"; then
		echo "  ln -s '$i'"
		ln -s "$i"
	fi
done
