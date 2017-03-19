#! /bin/sh
set -e
WHEEL_TOOL=`which wheel` /usr/bin/python2.7 setup.py sdist bdist_wheel
find dist -type f -exec gpg2 --detach-sign -a {} \;
twine upload dist/*
