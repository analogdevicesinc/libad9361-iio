#!/bin/bash -e

swdownloads_artifacts() {
    local linux_dist='Fedora-34 Fedora-28 Ubuntu-18.04 Ubuntu-20.04 Ubuntu-22.04 Debian-11 Debian-12 openSUSE-15.4 CentOS-7'
    for distribution in $linux_dist; do
        cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}/Linux-${distribution}"
        if [ "${distribution}" == "Fedora-34" ] || [ "${distribution}" == "Fedora-28" ] || [ "${distribution}" == "CentOS-7" ]; then
            find . -name '*.rpm' -exec mv {} ../"${distribution}_latest_${branch}_libad9361.rpm" ";"
        else
            find . -name '*.deb' -exec mv {} ../"${distribution}_latest_${branch}_libad9361.deb" ";"
        fi
        rm -r ../Linux-"${distribution}"
    done

	local macOS_dist='macOS-12 macOS-13-x64 macOS-13-arm64'
	for distribution in $macOS_dist; do
        cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}/${distribution}"
        find . -name '*.pkg' -exec mv {} ../"${distribution}_latest_${branch}_libad9361.pkg" ";"
        rm -r ../"${distribution}"
    done

	local windows_dist='2019 2022'
        for distribution in $windows_dist; do
		cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}"
                zip -r "Windows-VS-${distribution}-x64-latest_${branch}_libad9361".zip "Windows-VS-${distribution}-x64"
                rm -r "Windows-VS-${distribution}-x64"
        done

	local arm_dist='arm32v7 arm64v8 ppc64le x390x'
        for distribution in $arm_dist; do
                cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}/Ubuntu-${distribution}"
                find . -name '*.deb' -exec mv {} ../"Ubuntu-${distribution}_latest_${branch}_libad9361.deb" ";"
                rm -r ../Ubuntu-"${distribution}"
        done

	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}/Libad9361-Setup-Exe"
	mv libad9361-setup.exe ../libad9361-setup.exe
	rm -r ../Libad9361-Setup-Exe

}

branch=${2}
echo $branch
"${1}"_artifacts
