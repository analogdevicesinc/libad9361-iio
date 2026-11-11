
$COMPILER=$Env:COMPILER
$ARCH=$Env:ARCH
$BUILD_TYPE=if ($Env:CMAKE_BUILD_TYPE) { $Env:CMAKE_BUILD_TYPE } else { "Release" }

# Build type comes from CMAKE_BUILD_TYPE (passed by CI); defaults to Release.
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
	..

cmake --build . --config $BUILD_TYPE

if ( $LASTEXITCODE -ne 0 ) {
		throw "[*] cmake build failure"
	}

cp .\libad9361-iio.iss $env:BUILD_ARTIFACTSTAGINGDIRECTORY
