@rem for while this, for not taking time would be easier use distutils,
@rem but i dont like how it organize the folders anyway this...
@rem cd release
@rem del *
@rem cd ..
@rem del *.pyd
@rem make -f makefile
@rem now this one here works properly
copy ..\c\v.c .
swig.exe -python *.i
python setup.py build_ext --inplace