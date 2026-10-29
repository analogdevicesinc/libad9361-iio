#!/bin/bash -e

release_artifacts() {
	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}"

	# Bundle the Windows SDK artifacts (headers + libs) into a single zip.
	# There is normally one Windows-* directory (the VS build); handle any present.
	if ls -d Windows-*/ >/dev/null 2>&1; then
		mkdir -p Windows/include
		for d in Windows-*/; do
			d="${d%/}"
			# the libiio header is published once under include/
			if [ -f "${d}/iio.h" ]; then
				cp "${d}/iio.h" Windows/include/
				rm -f "${d}/iio.h"
			fi
			mv "${d}" Windows/
		done
		( cd Windows && zip -r Windows.zip ./* )
		mv Windows/Windows.zip ./
		rm -r Windows
	fi

	# Windows installer: flatten the setup .exe to the top level so it is
	# published (action-gh-release skips nested directories).
	if [ -d Libad9361-Setup-Exe ]; then
		find Libad9361-Setup-Exe -name '*.exe' -exec mv {} ./ ";"
		rm -r Libad9361-Setup-Exe
	fi

	# Flatten every remaining per-platform artifact directory to the top level.
	# Directories are discovered at runtime instead of hard-coded, so a matrix
	# change needs no edit here. Non-release bundles are left untouched.
	for d in */; do
		d="${d%/}"
		[ -d "${d}" ] || continue
		case "${d}" in
			Artifact-manifest|libad9361-dist|libad9361-libiio-upstream-*)
				continue
				;;
			macOS-*)
				# macOS .pkg names don't carry the arch, so append it.
				pkg=$(find "${d}" -name '*.pkg' | head -1)
				if [ -z "${pkg}" ]; then
					echo "ERROR: no .pkg artifact found in ${d}" >&2
					exit 1
				fi
				base=$(basename "${pkg}")
				mv "${pkg}" "${d}/${base%.pkg}-${d}.pkg"
				;;
		esac
		find "${d}" \( -name '*.deb' -o -name '*.rpm' -o -name '*.pkg' -o -name '*.tar.gz' \) -exec mv {} ./ ";"
		rm -r "${d}"
	done
}

check_artifacts() {
	cd "${BUILD_ARTIFACTSTAGINGDIRECTORY}"
	local manifest="artifact_manifest.txt"
	if [ ! -f "${manifest}" ]; then
		echo "ERROR: ${manifest} not found in ${BUILD_ARTIFACTSTAGINGDIRECTORY}"
		exit 1
	fi
	local missing=0
	while IFS= read -r line; do
		[ -z "${line}" ] && continue
		# Manifest paths may be recorded relative to the workspace
		# (artifacts/<tag>/<file>) or to the staging dir; accept either.
		rel="${line#artifacts/}"
		if [ -f "${rel}" ] || [ -f "${line}" ]; then
			echo "${line} exists."
		else
			echo "${line} does not exist."
			missing=$((missing + 1))
		fi
	done < "${manifest}"
	if [ "${missing}" -gt 0 ]; then
		echo "ERROR: ${missing} expected artifact(s) missing"
		exit 1
	fi
}

case "${1:-}" in
	release|check)
		"${1}"_artifacts
		;;
	*)
		echo "Usage: ${0} {release|check}" >&2
		exit 1
		;;
esac
