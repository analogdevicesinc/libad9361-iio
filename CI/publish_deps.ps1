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
	# The redist version folder (e.g. 14.4x.xxxxx) and the VC14x.CRT folder change
	# with the toolset shipped on the runner, so resolve the newest ones dynamically.
	$redistBase = 'C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Redist\MSVC'
	$version = Get-ChildItem -Path $redistBase -Directory |
		Where-Object { $_.Name -match '^\d+\.\d+' } |
		Sort-Object { [version]$_.Name } -Descending |
		Select-Object -First 1
	$crt = Get-ChildItem -Path (Join-Path $version.FullName 'x64') -Directory -Filter 'Microsoft.VC*.CRT' |
		Sort-Object Name -Descending |
		Select-Object -First 1
	cd $crt.FullName
	cp .\msvcp140.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
	cp .\vcruntime140.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
}

cd $src_dir

cp .\build\libxml2.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
cp .\build\libserialport-0.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
cp .\build\libusb-1.0.dll $env:BUILD_ARTIFACTSTAGINGDIRECTORY
