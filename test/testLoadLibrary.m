libName = 'libad9361';
hfile = '/usr/local/include/ad9361-wrapper.h';
loadlibraryArgs = {hfile,'includepath','/usr/local/include','addheader','ad9361.h'};
[a, b] = loadlibrary(libName, loadlibraryArgs{:});
libfunctions('libad9361')
unloadlibrary('libad9361')