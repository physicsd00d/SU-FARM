# ## ~~~~~~~~~~ Notes ~~~~~~~~~~~~~~~~ ##
# #
# #  Don't choose the name of your python class to be the same as a .cpp file you have already written
# #   because your source file will get overwritten by cython code for the python class!
# #
# #  extra_link_args=['']   # self explanatory
# #
# #  python setup.py build_ext --inplace
# #
# #
# #
# ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
# import os, time, copy

# import numpy # Do this so you can get the include directory

# from distutils.core import setup
# from distutils.extension import Extension
# from Cython.Distutils import build_ext

# from Cython.Build import cythonize

# # Contains cython wrapping information
# pyxFile = "CompactEnvelopeBuilder.pyx"

# # Names of all the source files that belong in the project
# #srcNames = ["Architecture.cpp", "Debris.cpp", "Footprint3D.cpp", "Point.cpp", "Random_Number.cpp", "timer.cpp", "Trajectory.cpp"]
# #srcNames = ["Architecture", "Debris", "Footprint3D", "Point", "Random_Number", "timer", "Trajectory"]
# srcNames = ["SkyGrid", "Debris", "Footprint3D", "Point", "Random_Number",  "timer", "PointCloud", "Trajectory"]
# #srcNames = ["Debris", "Footprint3D", "Point", "Random_Number"]

# # Location of the source files
# srcDir = os.path.relpath('cpp')

# # Anticipated location of object files
# objDir = os.path.realpath('../build/')

# # Generate paths to the files that can be passed to the compiler
# srcFiles = copy.copy(srcNames)
# objFiles = copy.copy(srcNames)
# compileThese = []
# linkThese = []

# for ix in range(len(srcNames)):
#     srcFiles[ix] = os.path.join(srcDir,srcNames[ix]) + ".cpp"
#     objFiles[ix] = os.path.join(objDir,srcNames[ix]) + ".o"

#     # Check which files have been modified and only compile them
#     if os.path.isfile(objFiles[ix]):    # Make sure the object file exists
#         objModTime = os.path.getmtime(objFiles[ix])
#         srcModTime = os.path.getmtime(srcFiles[ix])
        
# #        print time.ctime(objModTime) + " " + time.ctime(srcModTime) + "\n"
#         if (objModTime > srcModTime):
#             linkThese.append(objFiles[ix])  #This file has not been modified since it was last compiled, don't recompile it
#         else:
#             compileThese.append(srcFiles[ix])   #Timestamps don't match, this file was modified, recompile it
#     else:
#         compileThese.append(srcFiles[ix])   #The obj file doesn't exist, so it MUST be compiled


# # print srcFiles
# # print objFiles
# # print "\n\n\n"
# # print compileThese
# # print linkThese
# # print "\n\n\n"

# # Prepend the cython file
# compileThese.insert(0,pyxFile);





# extensions = [Extension(
#      "CompactEnvelopeBuilder",                 # name of extension
#      #                             ["CompactEnvelopeBuilder.pyx", "Trajectory.cpp", "timer.cpp", "Random_Number.cpp"], #  our Cython source
#      compileThese, #  our Cython source
#      include_dirs = [numpy.get_include(), '/sw/include/', srcDir],
#      libraries = ['gsl','gslcblas','kmlbase','kmldom','kmlengine'],
#      library_dirs = ['/sw/lib/'],
#     # points compiler to the correct std library (flags placed here OVERRIDE previous flags...e.g. -O0 here overrides default -O3)
#      extra_compile_args=['-mmacosx-version-min=10.8'],
#      extra_objects = linkThese,
#                 #export_symbols = ['CC=/Library/Frameworks/EPD64.framework/Versions/Current/bin/ccache-swig\ g++'],
#      language="c++"),]

# setup(
#     cmdclass = {'build_ext': build_ext},
#     name = "My hello app",
#     ext_modules = cythonize(extensions),
# )


# #
# # extensions = [
# #     Extension("CompactEnvelopeBuilder", sources = compileThese),
# #         # include_dirs = [...],
# #         # libraries = [...],
# #         # library_dirs = [...]),
# # ]
# # setup(
# #     # cmdclass = {'build_ext': build_ext},
# #     name = "My hello app",
# #     ext_modules = cythonize(extensions),
# # )

# # setup(
# #     # cmdclass = {'build_ext': build_ext},
# #     name = "My hello app",
# #     ext_modules = extensions,
# # )


# # setup(ext_modules= cythonize([Extension(
# #                              "CompactEnvelopeBuilder",                 # name of extension
# #                              #                             ["CompactEnvelopeBuilder.pyx", "Trajectory.cpp", "timer.cpp", "Random_Number.cpp"], #  our Cython source
# #                              compileThese, #  our Cython source
# #                              include_dirs = ['/sw/include/', srcDir],
# #                              libraries = ['gsl','gslcblas','kmlbase','kmldom','kmlengine'],
# #                              library_dirs = ['/sw/lib/'],
# #                             # points compiler to the correct std library (flags placed here OVERRIDE previous flags...e.g. -O0 here overrides default -O3)
# #                              extra_compile_args=['-mmacosx-version-min=10.8'],
# #                              extra_objects = linkThese,
# #                                         #export_symbols = ['CC=/Library/Frameworks/EPD64.framework/Versions/Current/bin/ccache-swig\ g++'],
# #                              language="c++")],  # causes Cython to create C++ source
# #       ))


# # setup(ext_modules= cythonize([Extension(
# #                              "CompactEnvelopeBuilder",                 # name of extension
# #                              #                             ["CompactEnvelopeBuilder.pyx", "Trajectory.cpp", "timer.cpp", "Random_Number.cpp"], #  our Cython source
# #                              compileThese, #  our Cython source
# #                              include_dirs = ['/sw/include/', srcDir],
# #                              libraries = ['gsl','gslcblas','kmlbase','kmldom','kmlengine'],
# #                              library_dirs = ['/sw/lib/'],
# #                             # points compiler to the correct std library (flags placed here OVERRIDE previous flags...e.g. -O0 here overrides default -O3)
# #                              extra_compile_args=['-mmacosx-version-min=10.8'],
# #                              extra_objects = linkThese,
# #                                         #export_symbols = ['CC=/Library/Frameworks/EPD64.framework/Versions/Current/bin/ccache-swig\ g++'],
# #                              language="c++")],  # causes Cython to create C++ source
# #       cmdclass={'build_ext': build_ext}))

# #setup(ext_modules=[Extension(
# #                   "CompactEnvelopeBuilder",                 # name of extension
# ##                             ["CompactEnvelopeBuilder.pyx", "Trajectory.cpp", "timer.cpp", "Random_Number.cpp"], #  our Cython source
# #                   compileThese, #  our Cython source
# #                   include_dirs = ['/sw/include/', srcDir],
# #                   libraries = ['gsl','gslcblas','kmlbase','kmldom','kmlengine'],
# #                   library_dirs = ['/sw/lib/'],
# #                   extra_compile_args=['-mmacosx-version-min=10.8'],     # points compiler to the correct std library
# #                    extra_objects = linkThese,
# ##                   export_symbols = ['CC=/Library/Frameworks/EPD64.framework/Versions/Current/bin/ccache-swig\ g++'],
# #                   language="c++")],  # causes Cython to create C++ source
# #      cmdclass={'build_ext': build_ext})

