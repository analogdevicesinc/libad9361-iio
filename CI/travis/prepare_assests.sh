#!/bin/bash -e

release_artifacts() {
	local linux_dist='Ubuntu-22.04 Ubuntu-24.04 Ubuntu-26.04 Debian-12 Debian-13 Fedora-42 Fedora-44 openSUSE-15.6 openSUSE-16.0'
	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}"
	for i in $linux_dist; do
		cd "Linux-${i}"
		find . -name '*.rpm' -exec mv {} ../ ";"
		find . -name '*.deb' -exec mv {} ../ ";"
		find . -name '*.tar.gz' -exec mv {} ../ ";"
		cd ../
		rm -r "Linux-${i}"
	done

	local macOS_dist='macOS-15-arm64 macOS-15-x64 macOS-26-arm64 macOS-26-x64 macOS-27-arm64'
	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}"
	for i in $macOS_dist; do
		cd "${i}"
		for pkg in *.pkg; do
			[ -f "$pkg" ] || continue
			base="${pkg%.pkg}"
			mv "$pkg" "${base}-${i}.pkg"
		done
		find . -name '*.pkg' -exec mv {} ../ ";"
		find . -name '*.tar.gz' -exec mv {} ../ ";"
		cd ../
		rm -r "${i}"
	done

	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}"
	mkdir -p Windows/include
	cp ./Windows-VS-2022-x64/ad9361.h ./Windows/include
	cd "Windows-VS-2022-x64"
	rm -f ad9361.h
	cd ../
	mv "Windows-VS-2022-x64" Windows/
	cd Windows
	zip -r ../Windows.zip ./*
	cd ../
	rm -r Windows

	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}/libad9361-Setup-Exe"
	find . -name '*.exe' -exec mv {} ../ ";"
	cd ../
	rm -r "libad9361-Setup-Exe"

	local arm_dist='Ubuntu-22.04-arm32v7 Ubuntu-22.04-arm64v8 Ubuntu-22.04-ppc64le Ubuntu-22.04-s390x Ubuntu-26.04-arm32v7 Ubuntu-26.04-arm64v8 Ubuntu-26.04-ppc64le Ubuntu-26.04-s390x Debian-12-arm64 Debian-12-armhf Debian-13-arm64 Debian-13-armhf'
	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}"
	for i in $arm_dist; do
		cd "${i}"
		find . -name '*.deb' -exec mv {} ../ ";"
		find . -name '*.tar.gz' -exec mv {} ../ ";"
		cd ../
		rm -r "${i}"
	done

	rm -rf "${BUILD_ARTIFACTSTAGINGDIRECTORY}/Artifact-manifest"
}

check_artifacts() {
	cd build
	while IFS= read -r line; do
		if [ -z "${line}" ]; then continue
		fi
		test -f ./artifacts/"${line}" && echo "${line} exist." || echo "${line} does not exist."
	done < "artifact_manifest.txt"
}

branch=${2}
echo $branch
"${1}"_artifacts
