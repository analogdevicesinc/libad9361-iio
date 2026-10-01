# https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_preference_variables?view=powershell-7.2#erroractionpreference
$ErrorActionPreference = "Stop"
$ErrorView = "NormalView"

$artifactsDir = "$env:BUILD_ARTIFACTSTAGINGDIRECTORY\Windows-VS-2022-x64"
dir $artifactsDir
# The .iss resolves its [Files] sources from ARTIFACTS_DIR / VCREDIST_DIR; both the
# libad9361/libiio DLLs and the VC runtime DLLs are staged into this one directory.
iscc "/DARTIFACTS_DIR=$artifactsDir" "/DVCREDIST_DIR=$artifactsDir" "$artifactsDir\libad9361-iio.iss"

Get-ChildItem $env:BUILD_ARTIFACTSTAGINGDIRECTORY -Force -Recurse | Remove-Item -Force -Recurse
cp C:\libad9361-setup.exe $env:BUILD_ARTIFACTSTAGINGDIRECTORY
