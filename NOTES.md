How to use cmake with pypi:
https://bloerg.net/2012/11/10/cmake-and-distutils.html

TODO
====

- Compile `_arv.so` as one unit, i.e. don't use shared library. It's just a
  mess with RPATH on OSX, and can better use LTO.
- Fix setup.py
- Use tox
