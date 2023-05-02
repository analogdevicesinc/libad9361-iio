
$COMPILER=$Env:COMPILER
$ARCH=$Env:ARCH

$src_dir=$pwd

choco install sphinx

if (!(Test-Path build)) {
	mkdir build
}

cp .\libad9361-iio.iss.cmakein .\build

cd build

cmake -G "$COMPILER" -A "$ARCH" `
        -DLIBIIO_LIBRARIES:FILEPATH=$pwd\libiio.lib `
        -DLIBIIO_INCLUDEDIR:PATH=$pwd `
        -DCMAKE_CONFIGURATION_TYPES=Release `
	..

cmake --build . --config Release

if ( $LASTEXITCODE -ne 0 ) {
		throw "[*] cmake build failure"
	}

cp .\libad9361-iio.iss $env:BUILD_ARTIFACTSTAGINGDIRECTORY
