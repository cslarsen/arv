cmake_minimum_required(VERSION 3.5)
project(arv-python CXX)

configure_file(arv/__init__.py arv/__init__.py COPYONLY)
configure_file(setup.py.in setup.py)
