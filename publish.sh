#! /bin/sh
set -e
WHEEL_TOOL=`which wheel` python setup.py sdist bdist_wheel
find dist -type f -exec gpg2 --detach-sign -a {} \;
twine upload dist/*
