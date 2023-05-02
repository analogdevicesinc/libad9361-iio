# https://docs.microsoft.com/en-us/powershell/module/microsoft.powershell.core/about/about_preference_variables?view=powershell-7.2#erroractionpreference
$ErrorActionPreference = "Stop"
$ErrorView = "NormalView"

dir $env:BUILD_ARTIFACTSTAGINGDIRECTORY\Windows-VS-2019-x64
iscc $env:BUILD_ARTIFACTSTAGINGDIRECTORY\Windows-VS-2019-x64\libad9361-iio.iss

Get-ChildItem $env:BUILD_ARTIFACTSTAGINGDIRECTORY -Force -Recurse | Remove-Item -Force -Recurse
cp C:\libad9361-setup.exe $env:BUILD_ARTIFACTSTAGINGDIRECTORY
