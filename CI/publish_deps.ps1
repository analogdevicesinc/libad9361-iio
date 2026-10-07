# https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_preference_variables?view=powershell-7.2#erroractionpreference
$ErrorActionPreference = "Stop"
$ErrorView = "NormalView"

$src_dir=$pwd
$COMPILER=$Env:COMPILER

if ($COMPILER -eq "Visual Studio 16 2019") {
	cd 'C:\Program Files (x86)\Microsoft Visual Studio\2019\Enterprise\VC\Redist\MSVC\14.29.30133\x64\Microsoft.VC142.CRT'
	cp .\msvcp140.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
	cp .\vcruntime140.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
}else {
	cd 'C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Redist\MSVC\14.29.30133\x64\Microsoft.VC142.CRT'
	cp .\msvcp140.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
	cp .\vcruntime140.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
}

cd $src_dir

cp .\build\libxml2.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
cp .\build\libserialport-0.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
cp .\build\libusb-1.0.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
