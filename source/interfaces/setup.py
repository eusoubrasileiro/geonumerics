#!/usr/bin/env python
"""
setup.py file for SWIG example
"""
from distutils.core import setup, Extension
from os import system

# put a call to swig from here before compiling
# create _wrap.c file
system('swig.exe -python *.i')
#since is very dificult to set the paths of source etc.. just put everything in
# the same folder and them compile
system('copy ..\c\dsprocessing.c .')
system('copy ..\c\dsprocessing.h .')

py_module = Extension('_dsprocessing', sources=['dsprocessing_wrap.c', 'dsprocessing.c'],)

setup (name = 'dsprocessing',
       version = '0.1',
       author      = "Andre L. Ferreira",
       description = """Simple swig example from docs""",
       ext_modules = [py_module],
       py_modules = ["dsprocessing"],
       )

# clean back what should not be in this folder
system('del dsprocessing.c dsprocessing.h')