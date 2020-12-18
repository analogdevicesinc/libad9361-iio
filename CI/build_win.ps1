
$COMPILER=$Env:COMPILER
$ARCH=$Env:ARCH

$src_dir=$pwd

if (!(Test-Path build)) {
	mkdir build
}
cd build

cmake -G "$COMPILER" -A "$ARCH" `
        -DLIBIIO_LIBRARIES:FILEPATH=$pwd\libiio.lib `
        -DLIBIIO_INCLUDEDIR:PATH=$pwd `
        -DCMAKE_CONFIGURATION_TYPES=Release `
	..

cmake --build . --config Release
