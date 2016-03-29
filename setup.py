## ~~~~~~~~~~ Notes ~~~~~~~~~~~~~~~~ ##
#
#  Don't choose the name of your python class to be the same as a .cpp file you have already written
#   because your source file will get overwritten by cython code for the python class!
#
#  extra_link_args=['']   # self explanatory
#
#  python setup.py build_ext --inplace
#
# Looks helpful
# http://stackoverflow.com/questions/33738885/python-setuptools-not-including-c-standard-library-headers
#
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
import sys, os, copy
import numpy # Do this so you can get the include directory
from setuptools import setup, Extension
from Cython.Build import cythonize

def readme():
    with open('../README.md') as f:
        return f.read()

# Don't accidentally install this package yet.  Only develop is allowed.
def catchArgs():
    args = sys.argv
    if args[1] != "develop":
        raise NotImplementedError("It is not advised that you do anything other than develop at this point.")

# Contains cython wrapping information
pyxFile = "CompactEnvelopeBuilder.pyx"
srcNames = ["Grid3D", "SkyGrid", "Debris", "Footprint3D", "Point", "Random_Number",  "timer", "PointCloud", "Trajectory"]

# Location of the Cython file (pyx)
pyxDir = os.path.realpath('src/') 
# Location of the C++ source files
cppDir = os.path.realpath('src/cpp') 
# Location of the C++ source files
pythonDir = os.path.realpath('src/python') 
# Anticipated location of object files (I'll put them here)
objDir = os.path.realpath('build/') 

# Generate paths to the files that can be passed to the compiler
srcFiles = copy.copy(srcNames)
objFiles = copy.copy(srcNames)
compileThese = []
linkThese = []

# Check which C++ files need to be compiled
for ix in range(len(srcNames)):
    srcFiles[ix] = os.path.join(cppDir,srcNames[ix]) + ".cpp"
    objFiles[ix] = os.path.join(objDir,srcNames[ix]) + ".o"

    # Check which files have been modified and only compile them
    if os.path.isfile(objFiles[ix]):    # Make sure the object file exists
        objModTime = os.path.getmtime(objFiles[ix])
        srcModTime = os.path.getmtime(srcFiles[ix])
        
#        print time.ctime(objModTime) + " " + time.ctime(srcModTime) + "\n"
        if (objModTime > srcModTime):
            linkThese.append(objFiles[ix])  #This file has not been modified since it was last compiled, don't recompile it
        else:
            compileThese.append(srcFiles[ix])   #Timestamps don't match, this file was modified, recompile it
    else:
        compileThese.append(srcFiles[ix])   #The obj file doesn't exist, so it MUST be compiled

# Prepend the cython file
compileThese.insert(0,pyxDir + '/' + pyxFile);

# Everything has been defined in terms of absolute paths, so I can change directories with no repurcussions
# Changing dir so that my .so and .egg will wind up here and the develop links will work
# os.chdir('build')

# List of Extension kwargs
# https://docs.python.org/2/distutils/apiref.html#distutils.core.Extension

compactEnvelopeModule = Extension(
     "CompactEnvelopeBuilder",                 # name of extension
     compileThese, #  our Cython source
     include_dirs = [numpy.get_include(), '/sw/include/', cppDir],
     libraries = ['gsl','gslcblas','kmlbase','kmldom','kmlengine'],
     library_dirs = ['/sw/lib/'],
    # points compiler to the correct std library (flags placed here OVERRIDE previous flags...e.g. -O0 here overrides default -O3)
     extra_compile_args=['-mmacosx-version-min=10.8'],
     extra_objects = linkThese,
                #export_symbols = ['CC=/Library/Frameworks/EPD64.framework/Versions/Current/bin/ccache-swig\ g++'],
     language="c++")


# Explains (kinda) how to call f2py from setuptools
# debrisPropModule = Extension(
#     "debrisPropagationTest",
#     ['moduleConstants.f90','modulePlanetEarth.f90','denEstimator.f90','linearInterpolation.f90','speedOfSound.f90','ode.o','debrisPropAllTime.f90'],
#     language="fortran")

# print "pythonDir = {0}".format(pythonDir)

# In the meanwhile, let's just say that you've gotta put the .so file into the python modules
import shutil
shutil.copy('build/DebrisProp/debrisPropagation.so', 'src/python/packages/FriscoLegacy/')
shutil.copy('build/DebrisProp/orbitProp.so', 'src/python/packages/FriscoLegacy/')

config = {
    'name': 'SU-FARM',
    # 'version': '0.1.0',
    # 'description': 'Air Pressure Packages',
    'author': 'Tom Colvin',
    'author_email' : 'tcolvin',
    'install_requires': ['scipy', 'numpy', 'pp'],
    'dependency_links': ['http://www.parallelpython.com/downloads/pp/pp-1.6.4.zip'],
    # 'packages': [pythonDir + '/Simulation'],
    'packages': ['Simulation','FriscoLegacy'],
    'package_dir': {'': 'src/python/packages'}, # Telling where to expect python packages
    # 'scripts': [],
    'ext_modules' : cythonize([compactEnvelopeModule])
    # 'long_description': readme()
}
# Try using dot notation on the packages?

catchArgs()  
setup(**config)











