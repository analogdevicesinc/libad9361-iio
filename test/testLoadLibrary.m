libName = 'libad9361';
hfile = '/usr/local/include/ad9361-wrapper.h';
loadlibraryArgs = {hfile,'includepath','/usr/local/include','addheader','ad9361.h'};
[a1, b1] = loadlibrary(libName, loadlibraryArgs{:});
libfunctions('libad9361')


hfile = '/usr/share/libiio/matlab/iio-wrapper.h';
loadlibraryArgs = {hfile,'includepath','/usr/local/include','addheader','iio.h'};
[a2, b2] = loadlibrary('libiio', loadlibraryArgs{:});
libfunctions('libiio')

% Create the network context
ip_address = '192.168.2.1';
iio_ctx = calllib('libiio', 'iio_create_network_context', ip_address);

% Check if the network context is valid
if (isNull(iio_ctx))
    iio_ctx = {}; %#ok<NASGU>
    unloadlibrary('libad9361')
    unloadlibrary('libiio')
    error('Could not connect to the IIO server!');
end

% Increase the object's instance count
fprintf('Connected to IP %s\n', ip_address);

nb_devices = calllib('libiio', 'iio_context_get_devices_count', iio_ctx);

% If no devices are present return with error
if(nb_devices == 0)
    unloadlibrary('libad9361')
    unloadlibrary('libiio')
    error('No devices were detected in the system!');
end
fprintf('Found %d devices in the system\n', nb_devices);

% Detect if the targeted device is installed
dev_found = 0;
for i = 0 : nb_devices - 1
    dev = calllib('libiio', 'iio_context_get_device', iio_ctx, i);
    name = calllib('libiio', 'iio_device_get_name', dev);
    fprintf('%s\n',name);
    if strcmp(name,'ad9361-phy')
        ret = calllib(libName,'ad9361_set_bb_rate',dev,int32(2e6));
        disp(ret);
    end
    clear dev;
end

iio_ctx = {};

unloadlibrary('libad9361')
unloadlibrary('libiio')