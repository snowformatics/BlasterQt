from distutils.core import setup
import py2exe
# cd D:\Users\lueck\Google Drive\Python\Software projects\BlasterQt
# python setup.py py2exe

setup(windows=[{"script": "blaster_main.py"}],

      options={"py2exe":
                {"includes": ['sip', 'decimal', 'PyQt4.QtSql', "PyQt4.QtCore", "PyQt4.QtGui", "PyQt4.QtNetwork"],
                  #"bundle_files":1,
                  #"optimize": 2,
                  "dll_excludes": ["mswsock.dll", "powrprof.dll"]}
                })

#options = {"build_exe" : {"includes" : "atexit" }}