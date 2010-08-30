#this makefile was made by hand trying to call the compiler
# avoiding to use distutils but its just a lose of time USE DISTIUTILS

# Tc is compile the especific file as a C file
flags=/c /nologo /Ox /MD /W3 /GS- /DNDEBUG
includes=-IC:\Python26\include -IC:\Python26\PC
outputdir=release
lflags=/DLL /nologo /INCREMENTAL:NO
libpaths=/LIBPATH:C:\Python26\libs /LIBPATH:C:\Python26\PCbuild
outputdirPyd=

#conclusion munifest not mandatory!!
all: swig dsprocessing_wrap dsprocessing link

swig:
	swig.exe -python dsprocessing.i

dsprocessing_wrap:
	cl.exe $(flags) $(includes) /Tcdsprocessing_wrap.c /Fo$(outputdir)\dsprocessing_wrap.obj

dsprocessing:
	cl.exe $(flags) $(includes) /Tcdsprocessing.c /Fo$(outputdir)\dsprocessing.obj

link:
	link.exe $(lflags) $(libpaths) /EXPORT:init_dsprocessing $(outputdir)\dsprocessing_wrap.obj $(outputdir)\dsprocessing.obj "/OUT:$(outputdirPyd)_dsprocessing.pyd" /IMPLIB:$(outputdir)\_dsprocessing.lib /MANIFESTFILE:$(outputdir)\_dsprocessing.pyd.manifest

manifest:
	mt.exe -nologo -manifest $(outputdir)\_dsprocessing.pyd.manifest -outputresource:$(outputdirPyd)_dsprocessing.pyd;2