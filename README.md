SU_FARM
=======

Coming Eventually...Stanford University Framework for Aircraft Risk Management

* I use a custom directory for system-wide packages that I install (/sw)
** Make sure that .profile / .bash_rc has /sw/bin prepended to the path and /sw/lib prepended to LD_LIBRARY_PATH
** Don't forget to actually export the variables after setting them, haha

* I had to install autoconf on my mac (it used to be included in previous versions of Xcode's tools)
** http://ftpmirror.gnu.org/autoconf/autoconf-2.68.tar.gz
** ./configure --prefix=/sw
** make
** sudo make install
* I had to do the same with automake
** http://ftp.gnu.org/gnu/automake/automake-1.15.tar.gz
* Cripes it never ends, I need the gnu libtool (even though osx comes with its own kind)
** http://ftpmirror.gnu.org/libtool/libtool-2.4.6.tar.gz
** If you still get errors about missing install-sh, try this: autoreconf -f -i -Wall,no-obsolete
* Must install https://code.google.com/p/libkml in order to be able to compile
** Actually the above code (in the downloads section) is really old and clang with throw errors.  If using clang, then get a newer version
** https://github.com/google/libkml
** But this version requires all the stuff I just had to install
** There are some unused variables that are causing the compilation to fail (-Werror).  I submitted a pull request to the repo so hopefully it'll get accepted.
** Had to manually copy the boost includes provided by libkml into a directory where they can be seen by the compiler: sudo cp -R boost /sw/include/
** If not install in one of the normal places, might need to add the proper include and lib paths in the cython setup.py


* Have to install GSL.  Now that I'm not using my own propagation code, the GSL is just for random numbers?

* TODO: Remove GSL dependency if possible
* TODO: Remove libkml and just write files in plain text
* TODO: Make the setup.py file actually place the object files in build.  Right now it compiles but the files dont go where i want them.
* Please excuse any curse words in the source code.

